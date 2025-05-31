# Creating a Python package
How to create a Python package for the bican-mccarroll-manuscript1 repository.
## Creating a package using PyCharm
1. Get started with GitHub by cloning the broadinstitute/bican-mccarroll-manuscript1 repository and creating a branch using the [basic GitHub workflow](https://docs.google.com/document/d/1A5T8ZVhNxP0GCQmmo1R8qaPPXah9cvSLI96U3TZdYos/edit?tab=t.0).
1. In PyCharm, go to `File` -> `New Project...`.
1. Set the location of the project to `<your-git-sandbox>/bican-mccarroll-manuscript1/python/bican_mccarroll_<your-package-name>`, where `<your-package-name>` is the name of your package.
1. I picked the following settings, which I'm not sure are the best, but pick a [Python version](https://devguide.python.org/versions/) that will be supported for at least two years:
   - Interpreter type: `Custom environment`
   - Environment: `Generate new`
   - Type: `Virtualenv`
   - Base python: `Python 3.12.x`
1. Click `Create`.
1. In the project view, create directories `src` and `tests` under the project root.
1. Right-click the `src` directory and select `New` -> `Python Package`. Name it `bican_mccarroll_<your_package_name>`.
1. Go to Settings -> Project -> Project Structure, and label 'src' and 'tests' as Sources and Tests, respectively.
1. Develop your package as usual, including unit tests as appropriate.
1. Put a `LICENSE` file in the root of your package, and a `README.md` file in the root of your package.  You can copy the ones in python/bican_mccarroll_helloworld/ and modify the README.md as appropriate.
1. Create a `pyproject.toml` file in the root of your package.  You can copy the one in python/bican_mccarroll_helloworld/ and modify it as appropriate.
1. When you are ready to commit your changes, use the [basic github workflow](https://docs.google.com/document/d/1A5T8ZVhNxP0GCQmmo1R8qaPPXah9cvSLI96U3TZdYos/edit?tab=t.0) to commit and push your changes.
