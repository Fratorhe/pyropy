=====================
Reader
=====================

Creating a new reader for your experimental facility is easy.
To do so, create a new Class that follows the propocol defined by :class:`pyropy.ExperimentReader`. 
In other words, create a class that requires as input the filename and folder, and it creates the objects time, temperature, rho and dRho. The way you create those is up to you, but they have to be created at instantiation (this is in the __init__ or in the __post_init__).

.. code-block:: python

    from pyropy import PyrolysisParallel

.. literalinclude:: ../../pyropy/tests/1React_parallel/data_parallel_verification.json
    :language: json