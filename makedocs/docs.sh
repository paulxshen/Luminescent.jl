cd "C:/Users/pxshe/OneDrive/Desktop/New folder"
chmod -R 775 .
rm -f docs/*
mv -f makedocs/build/* docs
 mv -f examples/*/*.mp4 docs/assets
