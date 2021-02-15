import casefoam
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

solutionDir = 'surfaces'
filename = 'p_freeSurf.raw'

def centreAndWall(caseComb, time, currentDataFrame):
    """Return the centre and wall pos.

    """
    centre = currentDataFrame.iloc[:, 1].min()
    wall = currentDataFrame.iloc[:, 1].max()
    df = pd.DataFrame(np.array([time, centre, wall], ndmin=2),
                      columns=['time','centre', 'wall'])
    print(df)
    df = df.set_index('time')

    return df

sol = casefoam.posField_to_timeSeries(solutionDir, filename, centreAndWall)
sol.sort_values('time',inplace=True)
sol[['centre','wall']].plot()
plt.show()
