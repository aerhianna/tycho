# Contributing

This document explains how to effectively write code to improve the TYCHO Batch Runner.

## Development Setup

### Ensure you have Python 3.8 or later installed.

```bash
$ python --version  # On your system the command might be `python3 --version`
Python 3.10.13
```

If you get an error, you will need to [install Python](https://www.python.org/downloads/).

### Run the unit test suite

If these tests pass, that means everything is working correctly out of the box.

```bash
make test
```

## Best Practices

For this project, it's strongly encouraged to stick with tools from the Python standard library to reduce the complexity of installing, using, and developing the batch runner.
  - e.g., using `argparse` instead of [`click`](https://click.palletsprojects.com/en/8.1.x/)
  - e.g. using `unittest` instead of [`pytest`](https://docs.pytest.org/en/8.0.x/)

Before implementing a new feature or fixing a bug, strive to write a unit test which exposes the shortcoming in the existing code. Then write the minimal amount of code required to address the issue, iterating as necessary until all tests pass. Repeat as necessary for each feature/bug.
  - This process, known in industry as "Test-Driven Development", helps ensure the software stays well tested and well organized (i.e. easy to modify later). It's not a foolproof way to build excellent software, but it helps a lot.