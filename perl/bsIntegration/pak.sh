#!/bin/sh

rm ~/bsIntegration.tgz
tar --one-file-system --exclude bsvir --exclude pak.sh --exclude "*.0*" --exclude ProteinogenicAA --exclude Casava --exclude bsI --exclude "t.*" --exclude "*.ini" -czvLf ~/bsIntegration.tgz .
ls -l ~/bsIntegration.tgz
