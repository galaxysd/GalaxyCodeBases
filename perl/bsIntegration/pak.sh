#!/bin/sh

tar --one-file-system --exclude bsvir --exclude pak.sh --exclude "*.0*" --exclude ProteinogenicAA --exclude Casava -czvLf ~/bsIntegration.tgz .
ls -l ~/bsIntegration.tgz
