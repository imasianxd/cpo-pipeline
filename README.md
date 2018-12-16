[![Build Status](https://travis-ci.org/Public-Health-Bioinformatics/cpo-pipeline.svg?branch=master)](https://travis-ci.org/Public-Health-Bioinformatics/cpo-pipeline)

# cpo-pipeline

An analysis pipeline for the purpose of investigating [Carbapenemase-Producing Organisms](https://en.wikipedia.org/wiki/Carbapenem-resistant_enterobacteriaceae).

This codebase is derived from the [CPO_Prediction](https://github.com/imasianxd/CPO_Prediction) workflow, written by [Justin Jia](https://github.com/imasianxd) in collaboration with [William Hsiao](https://github.com/wwhsiao) and [Matt Croxen](https://github.com/mcroxen).

## Requirements

It is expected that several bioinformatics tools are provided as tool-specific conda environments. Each environment is activated in turn by calling:

```
source activate <tool-version>
```

...and deactivated by calling:

```
source deactivate
```

In order for this to work, a conda installation must be available on the user's `PATH`.

If these environments are not already available on your system, create them as follows:

```
conda create -n <tool-version> tool=version
```

For example:

```
conda create -n mash-2.0 mash=2.0
```


