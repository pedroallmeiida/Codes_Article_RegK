library(reticulate)

### Instale essas bibliotecas antes de rodar
#py_install("pandas")
#py_install("mpmath")

## diretorio onde o seu python está instalado
use_python("C:\\Users\\pedro\\AppData\\Local\\r-miniconda\\envs\\r-python\\python.exe")

### Carregue o arquivo python
source_python("C:/Users/pedro/Dropbox/ARTIGO_REGRESSAO_K/Codes/meijer.py" , convert = TRUE )
