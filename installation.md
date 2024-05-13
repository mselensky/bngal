---
title: Installation
layout: default
has_children: false
permalink: /installation/
---

#### Table of contents
1. [Installing `bngal`](#installing-bngal)
2. [Command-line utility (`bngal-cli`)](#command-line-utility-recommended)
3. [R package only with CRAN](#r-package-only-with-cran)


# Installing `bngal`

Although `bngal` is released as a standalone R package and can be interactively used in an IDE such as RStudio, I strongly recommend running the [command-line utility](https://github.com/mselensky/bngal-cli) wrapper (`bngal-cli`) to simplify its use, especially for first-time users. `bngal-cli` is currently only tested on MacOS and Linux.

You can quickly install both the `bngal` R package and its command-line utility wrapper via the following instructions:

## Command line utility (recommended)

There are two ways to install the command line utility: DockerHub or Singularity. I recommend using one of the images hosted on [DockerHub](https://hub.docker.com/repository/docker/mjsel/bngal/tags?page=1&ordering=last_updated) (do note your chip architecture - discussed further below). An installation route is also available from [Anaconda](https://www.anaconda.com/products/distribution) or [CRAN](https://cran.r-project.org/).

### Docker (and Singularity)

Before pulling from DockerHub, please note your chip architecture - `arm64` (e.g., Apple Silicon) or `amd64 / x86_64` (e.g., Intel). You will want to pull the right image matched to your chip type. The `bngal` images are bootstrapped with [micromamba-docker](https://github.com/mamba-org/micromamba-docker):

Image Tag | Version | Architecture
:--- | :--- | :---
`mjsel/bngal:1.0.1` | `1.0.1` | `amd64 / x86_64`
`mjsel/bngal:1.0.1-arm64` | `1.0.1` | `arm64`
`mjsel/bngal:1.0.0` | `1.0.0` | `amd64 / x86_64`
`mjsel/bngal:1.0.0-arm64` | `1.0.0` | `arm64`

`bngal-cli` is easily installable if you use Docker. You will only need to install one of these, depending on your architecture:
```
docker pull mjsel/bngal:1.0.1
docker pull mjsel/bngal:1.0.1-arm64 
```

Alternatively, you can pull the same image if you use Singularity:

```
singularity pull docker://mjsel/bngal:1.0.1
singularity pull docker://mjsel/bngal:1.0.1-arm64
```

### Anaconda virtual environment
If you instead prefer to use conda, please follow the instructions as follows. 

1. [Install the appropriate Anaconda version](https://www.anaconda.com/products/distribution) for your operating system if you don't have it already.
2. Clone the `bngal-cli` GitHub repository into your directory of choice (`my-directory`) and run the setup script in a bash or zsh shell session. This will install the `bngal` R package within a new conda environment called "bngal":

{% highlight bash %}
cd my-directory
git clone https://github.com/mselensky/bngal-cli --branch v.1.0
cd bngal-cli
bash bngal-setup.sh
{% endhighlight %}

And that's it! Sit tight and grab a coffee while `bngal-cli` installs. It may take a few minutes.

Once you successfully install and activate the `bngal` environment, you can remove the `bngal-cli` folder. When the `bngal` environment is active, you will have access to two `bngal-cli` functions:

### `bngal-cli` functions

Function | Application
:--- | :---
`bngal-build-nets` | Build network model(s) according to defined cutoffs
`bngal-summarize-nets` | Summarize and visualize network statistics from `bngal-build-nets`

## R package only with CRAN

If you only want to use the `bngal` R package interactively, you can install it and its dependencies within an active R session via:
```
source("https://raw.githubusercontent.com/mselensky/bngal-cli/main/R/install-R-pkgs.R")
```
Please refer to the internal documentation when using the standalone R package.

