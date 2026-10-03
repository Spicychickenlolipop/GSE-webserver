import gzip
from pathlib import Path
import tempfile
import socket
import unittest
from unittest.mock import MagicMock, patch

import requests
import fetch_geo


class GeoDownloadTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.cache = Path(self.temp.name)
        self.target = self.cache / 'GSE12345_family.soft.gz'
        self.payload = gzip.compress(b'^SERIES = GSE12345\n!Series_title = test\n')

    def response(self, payload=None, status=200):
        response = MagicMock()
        response.__enter__.return_value = response
        response.status_code = status
        response.iter_content.return_value = [self.payload if payload is None else payload]
        session = MagicMock()
        session.__enter__.return_value = session
        session.get.return_value = response
        return session, response

    def test_https_download_then_local_parser(self):
        session, _ = self.response()
        with patch('fetch_geo.requests.Session', return_value=session), patch('fetch_geo.GEOparse.get_GEO') as parser:
            fetch_geo.fetch_gse(' gse12345 ', self.cache)
        self.assertEqual(session.get.call_args.args[0], 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE12nnn/GSE12345/soft/GSE12345_family.soft.gz')
        parser.assert_called_once_with(filepath=str(self.target), geotype='GSE', silent=True)
        self.assertEqual(self.target.read_bytes(), self.payload)
        self.assertEqual(list(self.cache.glob('*.part')), [])

    def test_valid_cache_never_downloads(self):
        self.target.write_bytes(self.payload)
        with patch('fetch_geo.requests.Session') as session, patch('fetch_geo.GEOparse.get_GEO'):
            fetch_geo.fetch_gse('GSE12345', self.cache)
        session.assert_not_called()

    def test_corrupt_cache_is_replaced(self):
        self.target.write_bytes(b'broken')
        session, _ = self.response()
        with patch('fetch_geo.requests.Session', return_value=session), patch('fetch_geo.GEOparse.get_GEO'):
            fetch_geo.fetch_gse('GSE12345', self.cache)
        self.assertEqual(self.target.read_bytes(), self.payload)

    def test_failed_download_is_not_cached(self):
        session, response = self.response()
        response.iter_content.side_effect = requests.ConnectionError('connection interrupted')
        with patch('fetch_geo.requests.Session', return_value=session), self.assertRaisesRegex(RuntimeError, 'over HTTPS'):
            fetch_geo.fetch_gse('GSE12345', self.cache)
        self.assertFalse(self.target.exists())
        self.assertEqual(list(self.cache.glob('*.part')), [])

    def test_404_has_actionable_error(self):
        session, _ = self.response(status=404)
        with patch('fetch_geo.requests.Session', return_value=session), self.assertRaisesRegex(ValueError, 'no public GEO SOFT file'):
            fetch_geo.fetch_gse('GSE12345', self.cache)
        self.assertFalse(self.target.exists())

    def test_wrong_series_is_rejected(self):
        session, _ = self.response(gzip.compress(b'^SERIES = GSE999\n'))
        with patch('fetch_geo.requests.Session', return_value=session), self.assertRaisesRegex(ValueError, 'does not contain series'):
            fetch_geo.fetch_gse('GSE12345', self.cache)
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_truncated_gzip_is_rejected(self):
        session, _ = self.response(self.payload[:-5])
        with patch('fetch_geo.requests.Session', return_value=session), self.assertRaisesRegex(RuntimeError, 'validate the download'):
            fetch_geo.fetch_gse('GSE12345', self.cache)
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_file_server_dns_failure_uses_main_geo_server(self):
        session, response = self.response(b'^SERIES = GSE12345\n!Series_title = test\n')
        failure = requests.ConnectionError('DNS failed')
        failure.__cause__ = socket.gaierror(11001, 'getaddrinfo failed')
        session.get.side_effect = [failure, response]
        with patch('fetch_geo.requests.Session', return_value=session), patch('fetch_geo.GEOparse.get_GEO'):
            fetch_geo.fetch_gse('GSE12345', self.cache)
        self.assertIn('www.ncbi.nlm.nih.gov/geo/query/acc.cgi', session.get.call_args.args[0])
        self.assertIn('targ=all', session.get.call_args.args[0])
        with gzip.open(self.target, 'rt') as handle:
            self.assertIn('^SERIES = GSE12345', handle.read())
        self.assertEqual(list(self.cache.glob('*.part')), [])

    def test_total_dns_outage_has_concise_actionable_message(self):
        session, _ = self.response()
        failure = requests.ConnectionError('internal pool details')
        failure.__cause__ = socket.gaierror(11001, 'getaddrinfo failed')
        session.get.side_effect = failure
        with patch('fetch_geo.requests.Session', return_value=session), self.assertRaisesRegex(RuntimeError, 'internet connection and DNS') as raised:
            fetch_geo.fetch_gse('GSE12345', self.cache)
        self.assertEqual(session.get.call_count, 2)
        self.assertNotIn('internal pool', str(raised.exception))
        self.assertEqual(list(self.cache.iterdir()), [])

    def test_fallback_html_page_is_not_cached_as_data(self):
        session, response = self.response(b'<html>Service unavailable</html>')
        session.get.side_effect = [requests.ConnectionError('primary offline'), response]
        with patch('fetch_geo.requests.Session', return_value=session), self.assertRaisesRegex(ValueError, 'does not contain series'):
            fetch_geo.fetch_gse('GSE12345', self.cache)
        self.assertEqual(list(self.cache.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
