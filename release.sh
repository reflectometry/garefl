#!/bin/sh

# Usage: ./release.sh
#
# Build an official release of garefl, tagged with the current date.
#
# To replace a previous version (e.g., from the previous day), use
#     VERSION=-yyyy.mm.dd ./release.sh
#

# Currently version is tied to date
VERSION="${VERSION:-`date +%Y.%m.%d`}"
echo "Creating garefl-$VERSION"
export VERSION

REFLWEB=reflectometry.org:web/danse
NCNRWEB=webster.ncnr.nist.gov:software
REPO=svn+ssh://svn@danse.us/reflectometry

function okay {
    echo $*
    echo -n "Press y to continue: "
    read ans
    test "$ans" != "y" && exit
}

# Check release notes
ls -l RELEASE-NOTES
head -10 RELEASE-NOTES
okay "Are the RELEASE-NOTES up to date, and tagged for $VERSION?"

svn update
svn status
okay "Are all files up to date in the repository?"

make clean
make check
okay "Was the build error free?"

make dist
ls release/*$VERSION*
okay "Are you ready to update the server?"
scp release/*$VERSION* $NCNRWEB/release
scp release/garefl.html $NCNRWEB
scp release/*$VERSION* $REFLWEB/download
svn copy $REPO/trunk/garefl $REPO/releases/garefl-$VERSION -m "tag release $VERSION"

head -20 RELEASE-NOTES
echo "Copy the release notes onto reflectometry.org and update the links"

