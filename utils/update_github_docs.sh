#!/bin/bash
# Usage: update_github_docs.sh ~/projects/loos-docs ~/project/LOOS "message goes here"
# THIS SCRIPT IS NOT THOROUGHLY TESTED, USE AT YOUR OWN RISK
# I'd suggest essentially doing this manually until we're sure the script is always correct

docs_dir=$1
loos_dir=$2
message=$3

cd $docs_dir
git pull
git checkout gh-pages
# This is unsafe as written, so BE VERY CAREFUL TO GET THE DIRECTORIES RIGHT
rm -rf *

cd $loos_dir/build
cmake --build . --target=docs

rsync -av html/* $docs_dir

cd $docs_dir

git add -A
git commit -m "$message [ci skip]"

git push
