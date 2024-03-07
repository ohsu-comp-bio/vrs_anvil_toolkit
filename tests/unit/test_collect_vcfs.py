import time

import pytest
from unittest.mock import patch, MagicMock
from vrs_anvil.collector import collect_manifest_urls


@pytest.fixture
def mock_manifest():
    manifest = MagicMock()
    manifest.vcf_files = ["http://example.com/file1.vcf", "gs://bucket/file2.vcf", "s3://bucket/file3.vcf", "file://my/local/file4.vcf"]
    manifest.num_threads = 20
    manifest.work_directory = "work/"
    return manifest


def download_google_blob_with_sleep(*args, **kwargs):
    time.sleep(1)
    return lambda bucket_name, source_blob_name, destination_file_name: "work/file2.vcf"


@patch('vrs_anvil.collector.download_http_file', return_value="work/file1.vcf")
@patch(target='vrs_anvil.collector.download_google_blob', new_callable=download_google_blob_with_sleep)
@patch('vrs_anvil.collector.download_s3_object', return_value="work/file3.vcf")
@patch('vrs_anvil.collector.create_symlink_to_work_directory', return_value="work/file4.vcf")
def test_collect_manifest_urls(
        mock_download_s3_object,
        mock_download_google_blob,
        mock_download_http_file,
        mock_create_symlink_to_work_directory,
        mock_manifest
):
    files = [_ for _ in collect_manifest_urls(mock_manifest)]
    print(files)
    # sort since the order of the files is not guaranteed
    assert sorted(files) == ['work/file1.vcf', 'work/file2.vcf', 'work/file3.vcf', 'work/file4.vcf']
