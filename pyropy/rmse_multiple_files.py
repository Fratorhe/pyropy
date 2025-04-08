import numpy as np
import spotpy


def rmse_multiple_files(evaluation, simulation):
    """
    Compute the Root Mean Squared Error (RMSE) for multiple files.

    The RMSE is calculated using the following formula:

    .. math::

     RMSE = \\sqrt{\\frac{1}{N}\\sum_{i=1}^{N}(e_{i}-s_{i})^2}

    Parameters
    ----------
    evaluation : list of lists
        Observed data to compare with simulation data.

    simulation : list of lists
        Simulation data to compare with evaluation data.

    Returns
    -------
    float
        Root Mean Squared Error (RMSE) value.
    """

    scale_coeff_dRho = 1000
    mses = []
    # Iterate through each objective and file pair
    for i, eval_data in enumerate(evaluation):
        for j, eval_file in enumerate(eval_data):
            mse_calc = spotpy.objectivefunctions.mse(eval_file, simulation[i][j])

            # Scale dRho (first objective) by the coefficient
            if i == 0:
                mse_calc *= scale_coeff_dRho

            mses.append(mse_calc)

    # Return the square root of the sum of the MSEs
    return np.sqrt(np.sum(mses))
