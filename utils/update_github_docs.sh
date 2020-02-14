#!/bin/bash
# Usage: update_github_docs.sh ~/projects/loos-docs ~/project/LOOS "message goes here"

docs_dir=$1
loos_dir=$2
message=$3

cd $docs_dir
git checkout gh-pages

cd $loos_dir
doxygen

rsync -av Docs/html $docs_dir

cd $docs_dir

git add -A
git commit -m "$message"

git push
