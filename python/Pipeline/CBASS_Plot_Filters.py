import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def PlotFilters (sFILTER):
#  Plot the filter for the enriched states
    inSizeGrid = np.ceil(np.sqrt(len(sFILTER))).astype(int)
    if inSizeGrid > 1:
        fig,ax3 = plt.subplots(inSizeGrid,inSizeGrid,figsize=(10,10))
        ax3 = ax3.ravel()
        count=0
        for idxFilters in range(inSizeGrid**2): #len(sFILTER)):
            if idxFilters >= len(sFILTER):
                ax3[count].set_axis_off()
            else: 
                ax3[count].imshow(sFILTER[idxFilters], aspect='auto')
                ax3[count].set_title('Filter {}'.format(count))
            count+=1
    else:
        fig,ax3 = plt.subplots(inSizeGrid,inSizeGrid,figsize=(5,5))
        for idxFilters in range(len(sFILTER)):
            ax3.imshow(sFILTER[idxFilters], aspect='auto')
            ax3.set_title('Filter')

    fig.tight_layout()
    plt.show()