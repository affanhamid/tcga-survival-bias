cat > install_gdc.sh << 'EOF'
#!/bin/bash
set -e

OS=$(uname -s)
ARCH=$(uname -m)

if [ "$OS" = "Darwin" ]; then
    if [ "$ARCH" = "arm64" ]; then
        URL="https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_OSX_x64-py3.8-macos-14.zip"
        INNER="gdc-client_2.3_OSX_x64.zip"
    else
        URL="https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_OSX_x64-py3.8-macos-12.zip"
        INNER="gdc-client_2.3_OSX_x64.zip"
    fi
elif [ "$OS" = "Linux" ]; then
    URL="https://gdc.cancer.gov/system/files/public/file/gdc-client_2.3_Ubuntu_x64-py3.8-ubuntu-20.04.zip"
    INNER="gdc-client_2.3_Ubuntu_x64.zip"
else
    echo "Unsupported OS: $OS"
    exit 1
fi

echo "Downloading gdc-client for $OS ($ARCH)..."
wget -q "$URL"
unzip -q "$(basename $URL)"
unzip -q "$INNER"
sudo mv gdc-client /usr/local/bin/
rm -f "$(basename $URL)" "$INNER"

echo "Done."
gdc-client --version
EOF

chmod +x install_gdc.sh