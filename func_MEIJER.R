library(reticulate)

### Install these libraries before running
#py_install("pandas")
#py_install("mpmath")

## directory where your python is installed
use_python("C:\\Users\\pedro\\AppData\\Local\\r-miniconda\\envs\\r-python\\python.exe")

### Load the python file
source_python("C:/Users/pedro/Dropbox/ARTIGO_REGRESSAO_K/Codes/meijer.py" , convert = TRUE )
