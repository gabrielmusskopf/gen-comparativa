import sys, os, traceback
# from PyQt5 import QtCore, QtGui, QtWidgets
# from PyQt5.QtWidgets import  *#QApplication, QMainWindow, QMessageBox, QRunnable
# from PyQt5.QtGui import *
# from PyQt5.QtWidgets import *
# from PyQt5.QtCore import *

from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg

from validate_email import validate_email

from ui_login import Ui_LoginWindow
from ui_methods import Ui_MethodWindow
from ui_loading import Ui_LoadingWindow
from ui_result import Ui_ResultWindow
from methods import *

from ui_splash_screen import Ui_SplashScreen


counter = 0


class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data
    
    error
        `tuple` (exctype, value, traceback.format_exc() )
    
    result
        `object` data returned from processing, anything

    progress
        `int` indicating % progress 

    '''
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(object)
    progress = pyqtSignal(int)



class Worker(QRunnable):
    '''	
    Classe que realiza o thread múltiplo, externo ao mesmo tempo que o resto do programa
    '''

    def __init__(self, fn, *args, **kwargs):
    	''' 
    	@param fn é a função de buca/alinhamento
    	'''
    	super(Worker,self).__init__()
    	self.fn = fn
    	self.args = args
    	self.kwargs = kwargs
    	self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        print("Pesquisa iniciada \n")

        try:
        	result = self.fn(*self.args, **self.kwargs)

        	if result == None:
        		traceback.print_exc()
	        	exctype, value = sys.exc_info()[:2]
	        	self.signals.error.emit((exctype, value, traceback.format_exc()))
	        
        except:		# Se der algum erro, vai printar no terminal o Traceback

        	print("-------------")
        	traceback.print_exc()
        	exctype, value = sys.exc_info()[:2]
        	self.signals.error.emit((exctype, value, traceback.format_exc()))
        	print("-------------\n\n")

        else:
        	self.signals.result.emit(result)	# Retorna o resultado

        finally:
        	self.signals.finished.emit()  # Done




class loginScreen(QMainWindow):
	'''
	Classe da janela de login, herda os métodos de QMainWindow
	'''
	
	def __init__(self):
		'''
		Setando a clsse com os objetos criados no arquivo .ui
		'''
		super().__init__()


		self.login_ui = Ui_LoginWindow()
		self.login_ui.setupUi(self)

		##
		# self.setWindowFlag(QtCore.Qt.FramelessWindowHint)
		# self.setAttribute(QtCore.Qt.WA_TranslucentBackground)


		# Instância da classe methodScreen
		self.mthdScrn = methodScreen(self.login_ui.emailEdit.text())
		# print(self.login_ui.emailEdit.text())'result_handle' referenced before assignment

		# Quando botão de login for clicado, conecta o SINAL ao SLOT
		self.login_ui.okButton.clicked.connect(self.checkEmail)

		self.setGIF()

		


	def setGIF(self):
		self.login_ui.movie = QMovie("images/dna02-unscreen.gif")
		self.login_ui.loginPicLabel.setMovie(self.login_ui.movie)
		self.login_ui.movie.setScaledSize(QSize(100,100))
		self.login_ui.movie.start()



	def checkEmail(self):  
	    self.isValid = CheckEmail(self.login_ui.emailEdit.text())
	    if self.isValid:
	    	self.loginToMethodScreen()



	def loginToMethodScreen(self):
		self.mthdScrn.show()
		self.close()
 



class methodScreen(QMainWindow):
	'''
	Classe da janela de métodos de pesquisa
	'''

	def __init__(self, email):
		'''
		Setando a clsse com os objetos criados no arquivo .ui
		'''
		super().__init__()
		self.email = email
		self.method_ui = Ui_MethodWindow()
		self.method_ui.setupUi(self)

		self.ldngScrn = loadingScreen(self)

		# Quando botão é clicado, vai para a função fileBrowser()
		self.method_ui.insertButton.clicked.connect(self.fileBrowser)
		# Quando botão de pesquisar, verifica se a pesquisa é válida
		self.method_ui.searchButton.clicked.connect(self.isValidSearch)



	def fileBrowser(self):
		'''
		Função definina no arquivo methods.py
		'''
		self.file_return = OpenDialogBox(self,self.method_ui) # Retorna só o ID	


	def isValidSearch(self):
		'''
		@param valid[0]: 0= inválido, 1 = válido
		@param valid[1]: 0 = local, 1 = web
		'''

		self.valid=IsValidSearch(self.method_ui)	# Retorna um objeto com [Pesquisa Válida, Método de pesquisa]

		if self.valid[0] == 1:
			self.methodToLoadingScreen()

			'''
			Se for válida, passa os param para a função arguments() da classe loading
			'''
			if self.valid[1] == 1:
				self.ldngScrn.arguments(1,self.method_ui,self.email)

			elif self.valid[1] == 0:
				self.ldngScrn.arguments(0,self.file_return,self.email)



	def methodToLoadingScreen(self):
		self.close()
		self.ldngScrn.show()




class loadingScreen(QMainWindow):
	'''
	Classe da janela de carregamento
	'''

	def __init__(self,methodObject):
		'''
		Configuração da janela com os objetos da .ui
		'''
		super().__init__()
		self.loading_ui = Ui_LoadingWindow()
		self.loading_ui.setupUi(self)

		self.mthdScrn = methodObject
		self.rsltScrn = resultScreen(methodObject)

		self.setGIF()

		'''
		Criando objeto thread para fazer a busca durante a tela de carregamento
		'''
		self.threadpool = QThreadPool()




	def setGIF(self):
		self.loading_ui.movie = QMovie("images/dna03-unscreen.gif")
		self.loading_ui.imageLabel.setMovie(self.loading_ui.movie)
		self.loading_ui.movie.setScaledSize(QSize(200,200))
		self.loading_ui.movie.start()



	def arguments(self, identifier, handler, email):
		'''
		Passando os param para a classe loading

		@param idenfifier: 0 = local, 1 = web
		@param handler: retorno da busca no Entrez
		'''
		self.identifier = identifier
		self.handler = handler
		self.email = email
		self.multiTrheadSearch()



	def multiTrheadSearch(self):
		'''
		Função para realizar a tarefa de busca sem interferir na aplicação normal do programa
		Realiza a função de alinhamento e só volta a interagir quando estiver completa (ou tenha erro)
		
		@param self.alignment função passada como parâmetro

		Após a thread ter terminado, os sinais são conectados as devidas funções
		'''
		worker = Worker(self.alignment)
		worker.signals.result.connect(self.shows)
		# worker.signals.progress.connect(self.progress)
		worker.signals.error.connect(self.error)
		worker.signals.finished.connect(self.thread_complete)
		self.threadpool.start(worker)



	def alignment(self):
		'''
		Função passada como param para o Thread

		@param identifier: 0 = local, 1 = web
		'''
		if self.identifier == 0: # Arquivo local
			resultAlignment = LocalAlignment(self.handler)

		elif self.identifier == 1:	# Arquivo web
			resultSearch = Search(self.handler,self.email) # Busca no banco de dados

			if resultSearch["IdList"] != []:
				resultAlignment = WebAlignment(resultSearch) # Realiza o alinhamento dessa busca

		return resultAlignment	# Retorna o resultado do alinhamento ou o erro ocorrido



	def shows(self, blast_record):
		'''
		Função conectada pelo sinal 'result' da thread
		'''

		if blast_record != None:
			'''
			Funções de exibições pertencentes a classe de resultados
			'''
			self.rsltScrn.showAlignments(blast_record)
			self.rsltScrn.showSites(blast_record)
			self.rsltScrn.showGraph(blast_record)
			self.loadingToResultWindow()



	def progress(self):
		'''
		Função conectada pelo sinal 'progress' da thread
		'''
		print("progress")



	def error(self, trbck):
		'''
		Função conectada pelo sinal 'error' da thread
		'''
		self.rsltScrn.returnToMethod()
		self.close()
		PopSearchError(trbck)


	def thread_complete(self):
		'''
		Função conectada pelo sinal 'finish' da thread
		'''
		print("THREAD COMPLETE!\n")


	def loadingToResultWindow(self):
		self.close()
		self.rsltScrn.show()




class resultScreen(QMainWindow):
	'''
	Classe da janela de resultados
	'''
	def __init__(self,methodObject):
		'''
		Passando os objetos da .ui
		'''
		super().__init__()

		self.mthdScrn = methodObject

		self.result_ui = Ui_ResultWindow()
		self.result_ui.setupUi(self)
		
		self.graph()

		self.result_ui.returnButton.clicked.connect(self.returnToMethod)



	def graph(self):
		self.result_ui.graphicsView = PlotWidget() #self.result_ui.basesFrame)
		self.result_ui.graphicsView.setMaximumSize(QtCore.QSize(850, 450))
		self.result_ui.graphicsView.setBackground('w')
		self.result_ui.graphicsView.setStyleSheet("QGraphicsView{\n"
"    color: rgb(0, 0, 0);\n"
"    background-color: rgb(40, 40, 40);\n"
"    border-radius: 5px;\n"
"}")
		self.result_ui.graphicsView.setObjectName("graphicsView")
		# self.result_ui.verticalLayout.addWidget(self.result_ui.graphicsView)

		return self.result_ui.graphicsView



	'''
	@param self.result_ui: objeto janela de resultados
	@param blast_record: objeto retornado no resultado positivo da Thread 
	@param self.seq_count: número de sequências a serem exibidas, valor inserido pelo usuário
	'''
	def showAlignments(self, blast_record):
		self.seq_count = int(self.mthdScrn.method_ui.sequencesAskLineEdit_2.text())
		ShowAlignments(self.result_ui, blast_record, self.seq_count)


	def showSites(self,blast_record):
		self.seq_count = int(self.mthdScrn.method_ui.sequencesAskLineEdit_2.text())
		ShowSites(self.result_ui, blast_record, self.seq_count)


	def showGraph(self, blast_record):
		self.seq_count = int(self.mthdScrn.method_ui.sequencesAskLineEdit_2.text())
		ShowGraph(self, blast_record, self.seq_count)

	# Ainda não
	def showPhylo(self):
		ShowPhylo(self,self.result_ui)
	


	def returnToMethod(self):
		'''
		Zerando os campos de texto e voltando para a tela de método
		'''
		self.mthdScrn.method_ui.searchEdit.setText('')
		self.mthdScrn.method_ui.insertEdit.setText('')
		self.mthdScrn.method_ui.sequencesAskLineEdit_2.setText('1')
		self.mthdScrn.show()
		self.close()
		


'''
Aplicação:
'''
if __name__ == "__main__":

	app = QtWidgets.QApplication(sys.argv)
	ui = loginScreen()	
	ui.show()
	sys.exit(app.exec_())