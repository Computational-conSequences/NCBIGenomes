curl -o /opt/biotools/NCBI2026/datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'
curl -o /opt/biotools/NCBI2026/dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/dataformat'
chmod +x /opt/biotools/NCBI2026/data*
xattr -dr com.apple.quarantine /opt/biotools/NCBI2026/datasets
xattr -dr com.apple.quarantine /opt/biotools/NCBI2026/dataformat
