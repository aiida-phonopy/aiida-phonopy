---
substitutions:
  AiiDA engine paper: '*Workflows in AiiDA: Engineering a high-throughput, event-based
    engine for robust and modular computational workflows*'
  AiiDA main paper: '*AiiDA 1.0, a scalable computational infrastructure for automated
    reproducible workflows and data provenance*'
  ESG: "Excellence Strategy of Germany\u2019s federal and state governments"
  Phonopy: '`phonopy`'
  Phonopy documentation: '`Phonopy` documentation'
  Phonopy paper: '*First principles phonon calculations in materials science*,'
  README.md of the repository: '`README.md` of the repository'
  aiida-core documentation: '`aiida-core` documentation'
  aiida-phonopy: '`aiida-phonopy`'
  aiida-quantumespresso: '`aiida-quantumespresso`'
  mapex: |-
    ```{image} images/MAPEX.jpg
    :width: 100%
    ```
  ubremen: |-
    ```{image} images/UBREMEN.png
    :width: 100%
    ```
---

```{eval-rst}
.. grid::
   :reverse:
   :gutter: 2 3 3 3
   :margin: 1 2 1 2

   .. grid-item::
      :columns: 12 4 4 4

      .. image:: images/logo_aiida.svg
         :width: 200px
         :class: sd-m-auto

   .. grid-item::
      :columns: 12 8 8 8
      :child-align: justify
      :class: sd-fs-5

      .. rubric:: AiiDA Phonopy

      An AiiDA plugin package to integrate the `Phonopy`_ software.
      Compute and store phonon related properties of materials with the popular open source `Phonopy`_ code
      with automatic data provenance provided by AiiDA.

      **aiida-phonopy version:** |release|

```

______________________________________________________________________

```{eval-rst}
.. grid:: 1 2 2 2
    :gutter: 3

    .. grid-item-card:: :fa:`rocket;mr-1` Get started
        :text-align: center
        :shadow: md

        Instructions to install, configure and setup the plugin package.

        +++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: installation/index
            :ref-type: doc
            :click-parent:
            :expand:
            :color: primary
            :outline:

            To the installation guides

    .. grid-item-card:: :fa:`info-circle;mr-1` Tutorials
        :text-align: center
        :shadow: md

        Easy examples to take the first steps with the plugin package.

        +++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: tutorials/index
            :ref-type: doc
            :click-parent:
            :expand:
            :color: primary
            :outline:

            To the tutorials

    .. grid-item-card:: :fa:`question-circle;mr-1` How-to guides
        :text-align: center
        :shadow: md

        Hands-on guides to achieve specific goals.

        +++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: howto/index
            :ref-type: doc
            :click-parent:
            :expand:
            :color: primary
            :outline:

            To the how-to guides

    .. grid-item-card:: :fa:`bookmark;mr-1` Topic guides
        :text-align: center
        :shadow: md

        Detailed background information on various concepts.

        +++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: topics/index
            :ref-type: doc
            :click-parent:
            :expand:
            :color: primary
            :outline:

            To the topic guides

    .. grid-item-card:: :fa:`cogs;mr-1` Reference guides
        :text-align: center
        :shadow: md

        Detailed reference guides on the application programming and command line interfaces.

        +++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: reference/index
            :ref-type: doc
            :click-parent:
            :expand:
            :color: primary
            :outline:

            To the reference guides
```

```{toctree}
:hidden: true
:maxdepth: 2

installation/index
tutorials/index
howto/index
topics/index
reference/index
```

# How to cite

If you use this plugin for your research, please cite the following works:

```{eval-rst}
.. highlights:: Atsushi Togo and Isao Tanaka, |Phonopy paper|_, Scripta Materialia **108**, 1-5 (2015)
```

```{eval-rst}
.. highlights:: Sebastiaan. P. Huber, Spyros Zoupanos, Martin Uhrin, Leopold Talirz, Leonid Kahle, Rico Häuselmann, Dominik Gresch, Tiziano Müller, Aliaksandr V. Yakutovich, Casper W. Andersen, Francisco F. Ramirez, Carl S. Adorf, Fernando Gargiulo, Snehal Kumbhar, Elsa Passaro, Conrad Johnston, Andrius Merkys, Andrea Cepellotti, Nicolas Mounet, Nicola Marzari, Boris Kozinsky, and Giovanni Pizzi, |AiiDA main paper|_, Scientific Data **7**, 300 (2020)
```

```{eval-rst}
.. highlights:: Martin Uhrin, Sebastiaan. P. Huber, Jusong Yu, Nicola Marzari, and Giovanni Pizzi, |AiiDA engine paper|_, Computational Materials Science **187**, 110086 (2021)

```

# Acknowledgements

We acknowledge support from:

```{eval-rst}
.. list-table::
    :widths: 60 40
    :class: logo-table
    :header-rows: 0

    * - The `U Bremen Excellence Chairs`_ program funded within the scope of the |ESG|_.
      - |ubremen|
    * - The `MAPEX`_ Center for Materials and Processes.
      - |mapex|
```

[aiida]: http://aiida.net
[aiida engine paper]: https://doi.org/10.1016/j.commatsci.2020.110086
[aiida main paper]: https://doi.org/10.1038/s41597-020-00638-4
[aiida-core documentation]: https://aiida.readthedocs.io/projects/aiida-core/en/latest/intro/get_started.html
[aiida-phonopy]: https://github.com/aiida-phonopy/aiida-phonopy
[aiida-quantumespresso]: https://github.com/aiidateam/aiida-quantumespresso
[esg]: https://www.dfg.de/en/research_funding/excellence_strategy/index.html
[mapex]: https://www.uni-bremen.de/en/mapex
[phonopy]: https://github.com/phonopy/phonopy
[phonopy documentation]: https://phonopy.github.io/phonopy/install.html
[phonopy paper]: http://dx.doi.org/10.1016/j.scriptamat.2015.07.021
[readme.md of the repository]: https://github.com/aiida-phonopy/aiida-phonopy/blob/develop/README.md
[u bremen excellence chairs]: https://www.uni-bremen.de/u-bremen-excellence-chairs
