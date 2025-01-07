# Tycho Batch Runner

Configure multiple Tycho simulations to run one after the other. Save time on resets/reconfigurations between runs! 

## Usage

### Ensure you have Python 3.8 or higher

```bash
$ python --version  # On your system the command might be `python3 --version`
Python 3.10.13
```

If you get an error, you will need to [install Python](https://www.python.org/downloads/).

### Install GhostScript

```bash
$ sudo apt install ghostscript
```

GhostScript converts PGPLOT's postscript output to PDFs. 

### Create a CSV file with one line per run you want to execute.

```csv
base_model,genex.fm,genex.fr,genex.omegbar,params.prefix,hrt.device,hrt.hrfile,metadata.label
a000000,0.7,0.7,4.34d-9,'D1','/ps','hr.D1',Simulation 1
a000000,0.8,0.8,5.42d-9,'D2','/ps','hr.D2',Simulation 2
a000000,0.7,1.0,5.75d-9,'D3','/ps','hr.D3',Simulation 3
```

`base_model` names the model file to use as the starting point for the simulation.

`genex.*` refers to the variables you want to set in the `genex.in` config file.
 - For example, `genex.zscale` will set the value of the variable `zscale` within `genex.in`

`params.*` refers to the variables you want to set in the `params.d` config file.
 - For example, `params.igraf` will set the value of the variable `igraf` within `params.d`

`hrt.*` refers to the variables you want to set in the `hrt.in` config file.
 - For example, `hrt.iyvar` will set the value of the variable `iyvar` within `hrt.in`

`metadata.*` are custom columns with data relevant to your organization. They are ignored by the batch runner.

> [!TIP]
>
> The CSV above is merely a sample CSV. You can use the batch CSV to update any field in the supported configs. You can also add additional metadata columns `metadata.*`

> [!NOTE]
> 
> No other config files are supported yet, but can be added easily to `ConfigUpdater.update`.

### Run the batch simulation

```bash
$ python run_batch.py --csvpath batch.csv --modeldir /path/to/models --rundir /path/to/rundir
```

This will...
 1. Open the csv at `csvpath`
 2. For each line in that CSV it will
    1. Copy the `base_model` from `modeldir` to `rundir/old.model`
    2. Update the configs in `rundir`
    3. Run `genex`
    4. Copy the `new.model` to `imodel`
    5. Run `tycho8`.

## Contributing

If the TYCHO Batch Runner is missing a feature or has a bug, please refer to the [CONTRIBUTING.md](./CONTRIBUTING.md) guide to learn the best practices for adding new code to this project.

## License

Copyright (c) 2024 - Joseph Hale. All Rights Reserved.

```
The TYCHO Batch Runner by Joseph Hale is licensed under the terms of the Mozilla
Public License, v 2.0, which are available at https://mozilla.org/MPL/2.0/.

You can download the source code for the TYCHO Batch Runner for free from
https://github.com/thehale/TYCHO.
```
<details>

<summary><b>What does the MPL-2.0 license allow/require?</b></summary>

### TL;DR

You can use files from this project in both open source and proprietary
applications, provided you include the above attribution. However, if
you modify any code in this project, or copy blocks of it into your own
code, you must publicly share the resulting files (note, not your whole
program) under the MPL-2.0. The best way to do this is via a Pull
Request back into this project.

If you have any other questions, you may also find Mozilla's [official
FAQ](https://www.mozilla.org/en-US/MPL/2.0/FAQ/) for the MPL-2.0 license
insightful.

If you dislike this license, you can contact me about negotiating a paid
contract with different terms.

**Disclaimer:** This TL;DR is just a summary. All legal questions
regarding usage of this project must be handled according to the
official terms specified in the `LICENSE` file.

### Why the MPL-2.0 license?

I believe that an open-source software license should ensure that code
can be used everywhere.

Strict copyleft licenses, like the GPL family of licenses, fail to
fulfill that vision because they only permit code to be used in other
GPL-licensed projects. Permissive licenses, like the MIT and Apache
licenses, allow code to be used everywhere but fail to prevent
proprietary or GPL-licensed projects from limiting access to any
improvements they make.

In contrast, the MPL-2.0 license allows code to be used in any software
project, while ensuring that any improvements remain available for
everyone.

</details>
