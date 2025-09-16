#!/bin/bash

# *****************************COPYRIGHT*******************************
# (C) Crown copyright 2024 Met Office. All rights reserved.
# For further details please refer to the file COPYRIGHT.txt
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************

github_url=${github_url:-"https://github.com/MetOffice/SimSys_Scripts.git"}
github_branch=${github_branch:-"main"}

echo "Git-scripts updater has started running"

rm -rf "${CYLC_WORKFLOW_SHARE_DIR}/imported_github_scripts"
if [[ $? != 0 ]]; then
    echo "Couldn't remove specified folder. Try checking permissions"
    exit 1
  else
    echo "Successfully removed old SimSys_Scripts git directory"
fi

git clone -b "$github_branch" --single-branch "$github_url" \
  "${CYLC_WORKFLOW_SHARE_DIR}/imported_github_scripts" 2>&1
if [[ $? != 0 ]]; then
    echo "Unable to clone remote git repo into specified location."
    echo "Check git branch, git url, destination path, permissions and git access"
    exit 1
  else
    echo "Git repo successfully cloned"
fi

echo "Copying suite_report and fcm_bdiff to cylc-run/bin"
cp "${CYLC_WORKFLOW_SHARE_DIR}/imported_github_scripts/suite_report.py" "${CYLC_WORKFLOW_RUN_DIR}/bin"
cp "${CYLC_WORKFLOW_SHARE_DIR}/imported_github_scripts/fcm_bdiff.py" "${CYLC_WORKFLOW_RUN_DIR}/bin"

echo "Github scripts updated"
