import numpy as np
import spotpy


def rmse_multiple_files(evaluation, simulation):
    """
    Root Mean Squared Error for more than 1 file

        .. math::

         RMSE=\\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2}

    :evaluation: Observed data to compared with simulation data.
    :type: list of lists

    :simulation: simulation data to compared with evaluation data
    :type: list of lists

    :return: Root Mean Squared Error
    :rtype: float
    """

    scale_coeff_dRho = 10000
    mses = []
    n_objectives = len(evaluation)  # 1 if 1 objective (dRho), 2 if 2 objectives (dRho, Rho)
    # n_objectives = 1  # 1 if only dRho
    n_files = len(evaluation[0])  # number of files to optimize
    for i in range(n_objectives):
        for j in range(n_files):
            mse_calc = spotpy.objectivefunctions.mse(evaluation[i][j], simulation[i][j])
            if i == 0:  # this means we are in dRho
                mses.append(mse_calc * scale_coeff_dRho)
            else:
                mses.append(mse_calc)
    return np.sqrt(sum(mses))