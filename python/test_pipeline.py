import gzip
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch, MagicMock

import numpy as np
import pandas as pd

from matrix_io import load_matrix, load_geo_matrix, validate_matrix
from run_pipeline import load_labels_from_file, infer_group_labels
import ml_models
import pipeline_methylation
import pipeline_rnaseq
import pipeline_microarray
from analysis_utils import welch_test
from fetch_geo import detect_data_type


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write(self, name, text):
        path = self.root / name
        if name.endswith('.gz'):
            with gzip.open(path, 'wt') as handle:
                handle.write(text)
        else:
            path.write_text(text, encoding='utf-8')
        return path

    def test_delimiters_and_compression(self):
        for sep in (',', '\t', ';', ' '):
            for ext in ('.txt', '.csv.gz'):
                with self.subTest(sep=sep, ext=ext):
                    path = self.write('matrix' + ext, sep.join(['gene', 's1', 's2']) + '\n' + sep.join(['g1', '1', '2']) + '\n')
                    self.assertEqual(load_matrix(path).shape, (1, 2))

    def test_series_matrix_comments(self):
        path = self.write('series.txt.gz', '!Series_title = test\n!series_matrix_table_begin\n"ID_REF"\t"GSM1"\t"GSM2"\n"a"\t1\t2\n!series_matrix_table_end\n')
        self.assertEqual(load_matrix(path).loc['a', 'GSM2'], 2)

    def test_invalid_inputs(self):
        for text in ('gene,s,s\na,1,2\n', 'gene,s\na,inf\n', 'gene,s\na,nope\n', 'gene,s\na,NA\n', 'gene,s\na,1\na,2\n'):
            with self.subTest(text=text), self.assertRaises(ValueError):
                load_matrix(self.write('bad.csv', text))

    def test_missing_features(self):
        matrix = load_matrix(self.write('missing.csv', 'gene,s1,s2\na,NA,NA\nb,1,NA\nc,2,3\n'))
        self.assertEqual(matrix.shape, (2, 2))
        self.assertTrue(np.isfinite(ml_models.run_pca(matrix.T)['coords']).all())

    def test_labels(self):
        path = self.write('labels.txt', 'sample,group\ns1,case\ns2,control\n')
        labels = load_labels_from_file(path, ['s2', 's1', 's3'])
        self.assertEqual(labels.iloc[0], 'control')
        self.assertTrue(pd.isna(labels.iloc[2]))
        with self.assertRaisesRegex(ValueError, 'No non-empty labels'):
            load_labels_from_file(path, ['other'])

    def test_methylation_missing_and_bounds(self):
        frame = pd.DataFrame([[.1, .2, .3, .4, np.nan], [.5, .6, .7, .8, .9]])
        result = pipeline_methylation.preprocess(frame)
        self.assertTrue(np.isfinite(result.to_numpy()).all())
        with self.assertRaisesRegex(ValueError, 'between 0 and 1'):
            pipeline_methylation.preprocess(pd.DataFrame([[2]]))

    def test_low_counts_and_zero_library(self):
        result = pipeline_rnaseq.preprocess(pd.DataFrame([[1, 0], [2, 0]]))
        self.assertTrue(np.isfinite(result.to_numpy()).all())
        with self.assertRaisesRegex(ValueError, 'positive counts'):
            pipeline_rnaseq.preprocess(pd.DataFrame([[0, 0]]))

    def test_microarray_does_not_relog(self):
        frame = pd.DataFrame([[2., 2.], [4., 4.]])
        with patch.object(pipeline_microarray, 'R_AVAILABLE', False):
            pd.testing.assert_frame_equal(pipeline_microarray.preprocess(frame), frame.T)

    def test_single_and_constant_samples(self):
        for frame in (pd.DataFrame([[1.]]), pd.DataFrame(np.ones((8, 3)))):
            pca = ml_models.run_pca(frame)
            self.assertTrue(np.isfinite(pca['coords']).all())
            self.assertEqual(pca['explained_variance_ratio'], [0., 0.])
            self.assertEqual(ml_models.run_clustering(frame)['k'], 1)

    def test_differential_constant_feature_does_not_poison_fdr(self):
        frame = pd.DataFrame({'variable': [1, 2, 6, 8], 'constant': [1, 1, 1, 1]}, index=list('abcd'))
        result = welch_test(frame, pd.Series(['x', 'x', 'y', 'y'], index=frame.index))
        self.assertEqual(result['gene'].tolist(), ['variable'])
        self.assertTrue(result.p_adj.between(0, 1).all())

    def geo(self, tables=None, urls=None):
        tables = tables or [pd.DataFrame(), pd.DataFrame()]
        return SimpleNamespace(gsms={f'GSM{i+1}': SimpleNamespace(table=table, metadata={'platform_id': ['GPL1']}) for i, table in enumerate(tables)},
                               metadata={'supplementary_file': urls or [], 'type': ['Expression profiling by array']}, gpls={})

    def test_geo_alternate_measurement_column(self):
        gse = self.geo([pd.DataFrame({'ID_REF': ['a'], 'counts': [1]}), pd.DataFrame({'ID_REF': ['a'], 'counts': [2]})])
        self.assertEqual(load_geo_matrix(gse, 'rnaseq').loc['a', 'GSM2'], 2)
        self.assertEqual(detect_data_type(gse), 'microarray')

    def test_supplementary_download_and_transpose(self):
        for content in (b'gene\tGSM1\tGSM2\na\t1\t2\n', b'sample,a,b\nGSM1,1,3\nGSM2,2,4\n'):
            response = MagicMock()
            response.__enter__.return_value = response
            response.status_code = 200
            response.iter_content.return_value = [content]
            with patch('requests.get', return_value=response):
                matrix = load_geo_matrix(self.geo(urls=['ftp://ftp.ncbi.nlm.nih.gov/geo/series/counts.txt']), 'rnaseq')
            self.assertEqual(matrix.loc['a', 'GSM2'], 2)

    def test_ambiguous_supplements(self):
        response = MagicMock()
        response.__enter__.return_value = response
        response.status_code = 200
        response.iter_content.return_value = [b'gene,GSM1,GSM2\na,1,2\n']
        with patch('requests.get', return_value=response), self.assertRaisesRegex(ValueError, 'Multiple processed matrices'):
            load_geo_matrix(self.geo(urls=['https://ftp.ncbi.nlm.nih.gov/a.txt', 'https://ftp.ncbi.nlm.nih.gov/b.txt']), 'rnaseq')

    def test_geo_sample_titles_and_annotation(self):
        gse = self.geo(urls=['https://ftp.ncbi.nlm.nih.gov/a.csv.gz'])
        gse.gsms['GSM1'].metadata['title'] = ['sample one']
        gse.gsms['GSM2'].metadata['title'] = ['sample two']
        response = MagicMock()
        response.__enter__.return_value = response
        response.status_code = 200
        response.iter_content.return_value = [gzip.compress(b'gene,annotation,sample one,sample two\n001,description,1,2\n')]
        with patch('requests.get', return_value=response):
            matrix = load_geo_matrix(gse, 'rnaseq')
        self.assertEqual(matrix.loc['001', 'GSM2'], 2)

    def test_classifier_with_partial_labels_and_missing_measurements(self):
        frame = pd.DataFrame({'a': [1, 2, np.nan, 10, 11, 12, 20], 'b': [3, 2, 1, 9, 10, 11, 20]}, index=list('abcdefg'))
        labels = pd.Series(['case'] * 3 + ['control'] * 3, index=list('abcdef'))
        result = ml_models.train_classifier(frame, labels)
        self.assertTrue(result['trained'])
        self.assertTrue(np.isfinite(result['cv_accuracy_mean']))

    def test_multiplatform_and_raw_only(self):
        gse = self.geo()
        with self.assertRaisesRegex(ValueError, 'Raw FASTQ'):
            load_geo_matrix(gse, 'rnaseq')
        gse.gsms['GSM2'].metadata['platform_id'] = ['GPL2']
        with self.assertRaisesRegex(ValueError, 'multiple platforms'):
            load_geo_matrix(gse, 'rnaseq')

    def test_ambiguous_phenotypes(self):
        pheno = pd.DataFrame({'characteristics': ['disease: case; sex: M'] * 2 + ['disease: control; sex: F'] * 2})
        self.assertIsNone(infer_group_labels(pheno, pheno.index))

    def test_reports_end_to_end(self):
        for data_type, values, scale in [('microarray', '4,5', 'auto'), ('rnaseq', '1,2', 'auto'), ('methylation', '.2,.8', 'auto'), ('rnaseq', '-1,2', 'processed')]:
            with self.subTest(data_type=data_type, scale=scale):
                matrix = self.write('input.csv.gz', 'gene,s1,s2\na,' + values + '\n')
                report = self.root / 'report.html'
                process = subprocess.run([sys.executable, str(Path(__file__).parent / 'run_pipeline.py'), '--matrix-file', str(matrix), '--data-type', data_type, '--data-scale', scale, '--job-id', 'test', '--out', str(report)], capture_output=True, text=True, timeout=120)
                self.assertEqual(process.returncode, 0, process.stderr)
                html = report.read_text(encoding='utf-8')
                data = json.loads(re.search(r'const pca = (.*);', html)[1])
                self.assertEqual(len(data['x']), 2)
                self.assertIn('Analysis notes', html)


if __name__ == '__main__':
    unittest.main()
