# Contributing Guide

Thank you for considering contributing to our Python project! By contributing, you help make our project better for everyone. Before you get started, please take a moment to review the following guidelines.

## Getting Started

1. **Fork the Repository:** Start by forking our project repository to your GitHub account. This will create a copy of the project under your account.

2. **Clone the Repository:** Clone the forked repository to your local machine using the following command:

```bash
git clone https://github.com/ohsu-comp-bio/vrs-python-testing
python3.10 -m venv venv; source venv/bin/activate
bash scripts/setup.sh
pip install -r requirements-dev.txt
pip install -e .
```
**Install Dependencies**: Ensure you have the necessary dependencies installed.

### Making Changes


Create a Branch: Before making changes, create a new branch for your feature or bug fix:

```bash
git checkout -b feature-name
```

Write Code: Make your code changes, keeping the coding style and project conventions in mind.

Write Tests: If applicable, write tests for your code changes to ensure they work as expected.

```bash
pytest tests/  #  --cov=<PACKAGE_NAME>

```

Check for Style: Run any code formatting tools or linters to maintain a consistent code style.

* pre commit tests

A reasonable set of checks, including running unit tests prior to each commit.  You can run these tests on demand by:

```
$ pre-commit install

$ pre-commit run --all-files
debug statements (python)................................................Passed
check python ast.........................................................Passed
fix utf-8 byte order marker..............................................Passed
check json...........................................(no files to check)Skipped
detect private key.......................................................Passed
check yaml...............................................................Passed
check for added large files..............................................Passed
check that scripts with shebangs are executable..........................Passed
check for case conflicts.................................................Passed
fix end of files.........................................................Passed
trim trailing whitespace.................................................Passed
mixed line ending........................................................Passed
run our unit tests.......................................................Passed

```

Commit Changes: Commit your changes with a clear and concise commit message:

```bash
git commit -m "Add feature X" -m "Fixes #123"
```

Push Changes: Push your changes to your forked repository:

```bash
git push origin feature-name
```
### Opening a Pull Request
Create a Pull Request: Open a pull request on the original repository. Provide a clear title and description of your changes.

Review Process: Participate in discussions and address feedback. Make additional commits if necessary.

Code Review: The project maintainers will review your code. Be prepared to make further changes if needed.

Merge: Once approved, your pull request will be merged. Congratulations!

Code of Conduct
Please note that our project has a Code of Conduct. We expect all contributors to adhere to its guidelines to ensure a positive and inclusive community.

Thank you for contributing to our project! Your efforts are highly appreciated. If you have any questions or need assistance, feel free to reach out to us.
