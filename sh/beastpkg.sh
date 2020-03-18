#!/usr/bin/env bash

arr=(
	"https://github.com/BEAST2-Dev/BEASTLabs/releases/download/v1.9.0/BEASTlabs.addon.v1.9.2.zip"
	"https://github.com/BEAST2-Dev/bModelTest/releases/download/v1.2.0/bModelTest.addon.v1.2.1.zip"
	"https://github.com/BEAST2-Dev/model-selection/releases/download/v1.5.0/MODEL_SELECTION.addon.v1.5.2.zip"
	"https://github.com/BEAST2-Dev/nested-sampling/releases/download/v1.1.0/NS.addon.v1.1.0.zip"
)

for ele in "${arr[@]}"; do
	name="$(basename "$ele")"
	root="$BEAST_PACKAGE_PATH/${name/.*/}"
	dest="$root/$name"
	curl -LJ0 --create-dirs -o "$dest" "$ele" && unzip -o "$dest" -d "$root" && rm -f "$dest"
done
