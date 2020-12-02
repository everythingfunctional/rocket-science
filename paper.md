---
title: 'Rocket Science: An Exercise in Refactoring Legacy Fortran'
tags:
  - Fortran
  - refactoring
  - rocket
authors:
  - name: Damian Rouson
    orcid: 0000-0002-2344-868X
    affiliation: 1
  - name: Brad Richardson
    orcid: 0000-0002-3205-2169
    affiliation: 1
  - name: Brian Laubacher
    affiliation: 2
affiliations:
 - name: Sourcery Institute
   index: 1
 - name: Autoliv
   index: 2
date: 2 December 2020
bibliography: paper.bib
---

# Summary

The near ubiquity of object-oriented programming support in modern programming languages
and the continued utility of legacy, procedural Fortran in important application domains
combine to make object-oriented design a common destination of code modernization efforts.
As Fortran has evolved into a multi-paradigm language, however, the peer-reviewed literature
contains far fewer detailed considerations of Fortran's support for functional programming.
This paper explores the symbiotic relationship between the two by demonstrating a legacy
code refactoring strategy in which the aim of writing pure functions motivates the choice
of abstractions and the ability to encapsulate function results facilitates the writing
of pure functions.  We show how a Fortran feature typically used for aliasing can be
applied to ensure the immutability of state.  We apply the resulting code refactoring
strategy to the modernization of an open-source pedagogical tool: a rocket motor
simulation mini-application.  We describe (1) the original, procedural Fortran 90 program
that uses global data, (2) a ground-up reimplementation of the same algorithms using
object-oriented design patterns in Fortran 2018, and (3) an evolutionary refactoring
that uses a consistent methodology for deriving purely functional abstractions from the
Fortran 90 code.  We contrast the two modern designs and demonstrate an evolutionary
path that could naturally lead the same resulting design for the ground-up reimplementation
and the evolutionary refactoring.

# Statement of need



# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

# Acknowledgements



# References
