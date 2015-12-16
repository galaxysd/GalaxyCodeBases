#!/bin/sh

echo "Before run me, uninstall every Adobe things and then use AdobeCreativeCloudCleanerTool."
read -p "Press [Enter] key to start cleaning..."

sudo rm -frv /Library/Application\ Support/Adobe/*
sudo rm -frv /Library/Application\ Support/regid.1986-12.com.adobe/
sudo rm -frv /Library/Logs/Adobe/

sudo rm -frv ~/Library/Application\ Support/Adobe/*
rm -v ~/Library/Application\ Support/AdobeWLCMCache.dat
rm -frv ~/Library/Caches/Adobe/*
rm -frv ~/Library/Caches/com.adobe.*
rm -frv ~/Library/Caches/CSXS/
rm -frv ~/Library/Caches/TemporaryItems/Adobe/
rm -frv ~/Library/Logs/Adobe/
rm -frv ~/Library/Preferences/Adobe*
rm -frv ~/Library/Preferences/com.{a,A}dobe.*
sudo rm -frv ~/Library/Saved\ Application\ State/com.adobe.*
rm -frv ~/Library/Logs/Adobe*

ADOBEERECEIPTS=$(pkgutil --pkgs=com.adobe*)
for ARECEIPT in $ADOBEERECEIPTS; do
  sudo pkgutil --forget $ARECEIPT
done

mkdir -p /Applications/Adobe/Adobe\ Acrobat\ DC/
