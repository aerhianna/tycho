# Copyright (c) 2024 Joseph Hale, Laura Pang
# 
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from argparse import ArgumentParser
from typing import Dict, List, TypedDict
from csv import DictReader
from pathlib import Path
import shutil
import re
import subprocess
import os


class Args(TypedDict):
    csvpath: str
    rundir: str
    modeldir: str


def parse_args(args: List[str]) -> Args:
    parser = ArgumentParser()
    parser.add_argument("--csvpath", required=True)
    parser.add_argument("--rundir", required=True)
    parser.add_argument("--modeldir", required=True)
    namespace = parser.parse_args(args)
    return {
        "csvpath": namespace.csvpath,
        "rundir": namespace.rundir,
        "modeldir": namespace.modeldir,
    }


def parse_csv_line(line: dict):
    parsed = {}
    for key, value in line.items():
        prefix, name = key.split(".") if "." in key else ("", key)
        parsed[prefix] = {**parsed.get(prefix, {}), name: value}
    return parsed


class ConfigFile(dict):
    configtext: str
    _regexes: Dict[str, re.Pattern]

    def __init__(self, configtext: str):
        self.configtext = configtext
        self._regexes = {}

    def __getitem__(self, key: str) -> str:
        result = self.__read(key)
        if result:
            return result[2]  # line, var, value = result
        else:
            raise KeyError(key)

    def __setitem__(self, key: str, value: str):
        result = self.__read(key)
        if result:
            newtext = f"{result[1]}{value}"
            self.configtext = re.sub(self.__regex(key), newtext, self.configtext)
        else:
            raise KeyError(key)

    def __str__(self) -> str:
        return self.configtext

    def __regex(self, key: str):
        if key not in self._regexes:
            self._regexes[key] = re.compile(rf"('{key}'\s+)(.*)")
        return self._regexes[key]

    def __read(self, key: str):
        return self.__regex(key).search(self.configtext)


class ConfigUpdater:
    updates: dict
    rundir: str
    modeldir: str

    def __init__(self, updates: dict, rundir: str, modeldir: str):
        self.updates = updates
        self.rundir = rundir
        self.modeldir = modeldir

    def update(self):
        if "" in self.updates:
            self._update_root(self.updates[""])
        if "genex" in self.updates:
            self._update_genex(self.updates["genex"])
        if "params" in self.updates:
            self._update_params(self.updates["params"])
        if "hrt" in self.updates:
            self._update_hrt(self.updates["hrt"])

    def _update_root(self, updates: dict):
        if "base_model" in updates:
            src = Path(self.modeldir, updates["base_model"])
            dst = Path(self.rundir, "old.model")
            print(f"MODEL COPY: {src} -> {dst}")
            shutil.copy(src, dst)

    def _update_genex(self, updates: dict):
        self._update_config("genex.in", updates)

    def _update_params(self, updates: dict):
        self._update_config("params.d", updates)

    def _update_hrt(self, updates: dict):
        self._update_config("hrt.in", updates)

    def _update_config(self, config_name: str, updates: dict):
        config_path = Path(self.rundir, config_name)
        print(f"CONFIG UPDATE: {config_path}")

        with open(config_path) as f:
            config = ConfigFile(f.read())

        for key, value in updates.items():
            print(f"SET {key} = {value}")
            config[key] = value

        with open(config_path, "w") as f:
            f.write(str(config))


def run_simulation(rundir: str, updates: dict):
    os.chdir(rundir)
    subprocess.run(["./genex"])
    shutil.copy("new.model", "imodel")
    subprocess.run(["./tycho8"])
    subprocess.run(["./hrt"])
    subprocess.run(['ps2pdf','pgplot.ps', updates['params']['prefix']+'.pdf'])

if __name__ == "__main__":
    import sys

    args = parse_args(sys.argv[1:])
    with open(args["csvpath"], "r") as f:
        csv = DictReader(f)
        for line in csv:
            updates = parse_csv_line(line)
            ConfigUpdater(updates, args["rundir"], args["modeldir"]).update()
            run_simulation(args["rundir"], updates)
