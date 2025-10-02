How the documentation works
===========================

Structure of the documentation:

* ``./doc`` holds all documentation-related information.
    * ``build`` holds the pretty html version of the documents once sphinx builds
      them.
    * ``source`` holds the ``.rst`` markdown documentation files.
        * ``conf.py`` gives Sphinx settings to build the documentation.
        * ``index.rst`` is the root document, which holds the table of contents
          of all rst files which will be in the built documentation (without the
          file extension). Any files not in this list will not be included, and
          Sphinx will warn you that they are not included when the documentation
          is built.
        * ``spacedoc.rst`` contains the documentation which is autobuilt from the
          docstrings in ``/SPaCE.py``, which should describe how the SPaCE package
          works.
        * ``tutorial.rst`` is essentially a clarified version of ``/SPaCE examples.ipynb``,
          reworded and formatted for Sphinx. Any updates to the examples should also be
          reflected here.
        * ``metadoc.rst`` is this page.
        * ``_static/`` is a directory for images which need to be included in the
          documentation, which can then be incorporated with an image directive:

.. code:: rst

    .. image:: _static/image.png

For a description of how to use the ReST format, see the [official documentation](https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html)
