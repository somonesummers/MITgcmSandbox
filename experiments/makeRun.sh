#!/usr/bin/env bash
set -e
TIME="$(date +"%H%M%S")" 
unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     MACHINE="Linux";;
    Darwin*)    MACHINE="Mac";;
esac
echo "Idenitfied machine as ${MACHINE}"
echo "Already build, clean up run folder, then make simlinks and run"
cd results
touch test.txt #this ensures the dir is not empty
rm *
ln -s ../input/* .
cp ../build/mitgcmuv .

echo "Ready to run -=三(ง ˙o˙)ว"
if [ "$MACHINE" == "Mac" ];
then
	time ./mitgcmuv > Report$TIME.txt
else
	./mitgcmuv
fi
echo "Done running  ദി(˵•̀ᴗ-˵)✧"
