## SYSTEM IMPORT
import os, sys, traceback

## PYQT5 IMPORTS
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import QMovie, QTransform
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox, QFileDialog
from PyQt5.QtCore import *

from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg 

## BIOPYTHON IMPORTS
import Bio
from Bio import Entrez, SeqIO, AlignIO, Phylo  
from Bio.Phylo import PhyloXMLIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq

## OTHER USEFUL IMPORTS
from validate_email import validate_email
from datetime import datetime
import math



# Checa se tem formato de email
def CheckEmail(email):  
    if email == 'admin':
        Entrez.email = 'gabrielgmusskopf@gmail.com'
        print("Bem-vindo(a), ", Entrez.email)
        return True
    if( validate_email( email )):  
        print(email)
        return True
    else:
    	print("nop")
    	NotValidEmailPopUp()



def NotValidEmailPopUp():
	pop = QMessageBox()
	pop.setWindowTitle("Erro")
	pop.setText("Por favor, insira um email válido")
	pop.setIcon(QMessageBox.Critical)
	pop.exec_()



def PopMultipleMethods(self):
    pop = QMessageBox()
    pop.setWindowTitle("Erro")
    pop.setText("Mais de um método selecionado")
    pop.setIcon(QMessageBox.Critical)
    pop.exec_()



def PopSearchError(error):
    pop = QMessageBox()
    pop.setWindowTitle("Erro")
    pop.setText("Erro na pesquisa ou inserção do arquivo\n" + str(error))
    pop.setIcon(QMessageBox.Critical)
    pop.exec_()



# Janela de pesquisa de arquivo
def OpenDialogBox(self,method_ui):
    options = QFileDialog.Options()
    filename = QFileDialog.getOpenFileName(filter="Fasta files (*.fasta)", options=options)
    path = filename[0]
    method_ui.insertEdit.setText(path)


    if path != '':
        record = SeqIO.read(path,"fasta")
        return record.id
    else:
        return "None"



def IsValidSearch(self):
    '''
    [0] = 1 se for válido
    [1] = 0 se for local
    [1] - 1 se for web 
    '''

    if (self.searchEdit.text() != "" or self.insertEdit.text() != ""):
        if (self.searchEdit.text() != "" and self.insertEdit.text() == ""): # Se for web
            return [1,1]
        elif (self.searchEdit.text() == "" and self.insertEdit.text() != ""): # Se for local
            return [1,0]
        elif (self.searchEdit.text() != "" and self.insertEdit.text() != ""): # Dois métodos selecionados
            PopMultipleMethods(self)
            return [0,0]
    else:
        PopSearchError(self)
        return [0,0]



def Search(method_ui,email):
    Entrez.email = email

    '''
    ~ Módulo Entrez
    
    Retorna um XML

    Exemplo:
    >>> from Bio import Entrez
    >>> Entrez.email = "Your.Name.Here@example.org"
    >>> handle = Entrez.einfo() # or esearch, efetch, ...
    >>> record = Entrez.read(handle)
    >>> handle.close()

    ~ Entrez.esearch()
    Busca e retorna os IDs

    Em caso de erro de rede, retorna IOError exception

    Parâmetros requeridos:
    @param db: database 

    Para mais parâmetros: https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch 
    '''

    term = method_ui.searchEdit.text()

    if ' ' in term:
        term = term.replace(' ','+')

    try:
        handle = Entrez.esearch(db="nucleotide", term = term, idtype="acc", retmax = 1) # Retorna um XML
        record = Entrez.read(handle) #lendo as infos geradas pela pesquisa # Converte XML para estrutura de dados python (dicionários)
        
        # print(record.keys())
        # print(record)

        h = Entrez.efetch(db="nucleotide", id=record['IdList'], rettype="fasta", retmode="text")
        record_fasta = SeqIO.read(h,'fasta')   #SeqIO.read() somente um ID
        # print(seq.id)
        # print(len(seq.seq))

        handle.close()
    except:
        return False

    return record_fasta
    


def LocalAlignment(id):

    print("Entrou no alinhamento de arquivo local: ",datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("ID de alinhamento: ", id)


    try:
        result_handle = NCBIWWW.qblast("blastn", "nt", id)
        blast_record = NCBIXML.read(result_handle) 

        result_handle.close()
    except:
        traceback.print_exc()
        return None


    # try:
    #     with open(path,"w") as save_file: # Para automaticamente fechar o "save_file"
    #         save_file.write(result_handle.read())
    #         result_handle.close()
    
    # except IOError:
    #     print("Erro na escrita do arquivo de alinhamento: " + path)


    #########
    # result_handle = open("/results/alignment_result_test.xml","r")
    #########

    # path = "/results/alignment_result.xml"

    # try:
    #     with open(path,"r") as result_handle:
    #         blast_record = NCBIXML.read(result_handle)

    # except:
    #     return None


    print("Fim do alinhamento: " + datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second + "\n")
    return blast_record



def WebAlignment(result_search):
    '''
    ~ Módulo Blast

    Função qblast() é responsável por chamar a versão online do Blast

    Parâmetros requeridos:

    @param 'blastn' programa blast para usar   >>>  https://blast.ncbi.nlm.nih.gov/Blast.cgi
    @param 'nt' database para buscar           >>>  ftp://ftp.ncbi.nlm.nih.gov/pub/factsheets/HowTo_BLASTGuide.pdf
    @param id pode ser o id da query ou também a sequência no formato fasta

    Para mais parâmetros:
    >>> from Bio.Blast import NCBIWWW
    >>> help(NCBIWWW.qblast)
    '''

    print("######################")
    print("Entrei no alinhamento da pesquisa web:", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("ID de alinhamento: ", result_search.id)

    # id = result_search.id
    # try:
    #     result_handle = NCBIWWW.qblast("blastn", "nt", id)
    #     blast_record = NCBIXML.read(result_handle) 
    #     result_handle.close()
    # except:
    #     traceback.print_exc()
    #     return None



    ########### Test ###########

    path = "results/alignment_result_test.xml"
    # path_test = "results/alignment_result_test.xml"

    try:
        with open(path,"r") as result_handle:
            # exec('blast_record = NCBIXML.read(result_handle)')
            blast_record = NCBIXML.read(result_handle)

    except:
        traceback.print_exc()
        return None

    #############################


    print("Fim alinhamento: ", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("######################\n")

    return blast_record


'''
Função que exibe alinhamentos na tela

@param self classe de onde foi chamado
@param results lista com o resultado da pesquisa lido por SeqIO e objeto em formato de dicionário que retornou do alinhamento
@param seq_count número de sequências que aparecerão na tela


Nas funções de exibição foi necessário escrever primeiro em um arquivo, para que possa 
posteriormente ser posto na label da janela criada

Exemplo de acesso à informações no blast_record:

>>> try:
>>>     for alignment in blast_record.alignments:
>>>         for hsp in alignment.hsps:
>>>             print("Score: ", hsp.score)
>>>             print("Bits: ", hsp.bits)
>>>             print("E: ", hsp.expect)
>>>             print("num_alignments: ", hsp.num_alignments)
>>>             print("Identities: ", hsp.identities)
>>>             print("Positives: ", hsp.positives)
>>>             print("Gaps: ", hsp.gaps)
>>>             print("Strand: ", hsp.strand)
>>>             print("Frame: ",hsp.frame)
>>>             print("Query start: ",hsp.query_start)
>>>             print("Sbjct star: ", hsp.sbjct_start)
>>> except:
>>>     traceback.print_exc()
>>>     pass

'''
def ShowAlignments(self, results, seq_count):


    print("###################")
    print("Entrou na função de exibir o alinhamento:", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)

    search_record = results[0]
    blast_record = results[1]

    query = search_record.seq

    count = 0
    maxs = 0
    queryies=[]
    matches=[]
    subjects=[]
    lengths = []


    ### Para escrever na tela, é preciso gravar em um arquivo ###

    path = "results/result.txt"

    try:
        with open(path ,'w') as result_file:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if count < seq_count:   
                        count+=1

                        ########### Escreve abeçalho ###########
                        result_file.write("\n *** Alinhamento *** \n" +
                        "Sequência: " + str(alignment.title) + "\n" + 
                        "Comprimento: "+ str(alignment.length) + "\n"+
                        "Número de bases diferentes: " + str(hamming(hsp.query,hsp.sbjct))+"\n\n")

                        ########### Escreve sequências ###########
                        # Para escrever do mesmo jeito que no BLAST
                        # math.ceil() arredonda o número para cima
                        for b in range( math.ceil( len( hsp.query ) / 60 )):


                            result_file.write(
                                str(hsp.query_start + (61 * b ) ) + " " + str(hsp.query[ b*60 : (b+1)*60 ] ) + " " + str((hsp.query_start + 60) + (b*61)) + "\n" + 
                                str(hsp.sbjct_start + (61 * b ) ) + " " + str(hsp.sbjct[ b*60 : (b+1)*60 ] ) + " " + str((hsp.sbjct_start + 60) + (b*61)) + "\n\n")

                        result_file.write("\n\n")


                        if alignment.length > maxs:
                            maxs = alignment.length
    except:
        traceback.print_exc()

    ''' Setando o tamanho conforme a quantidade de sequências '''
    self.alignmentScrollFrame.setMinimumSize(QtCore.QSize(1200,  25000)) #count*15000)) 


    try:
        with open(path,"r") as r:
            self.alignmentResultLabel.setText(r.read())

    except IOError:
        print("Erro na leitura: " + path)


    # print("\n" + "Tem ",count,"sequências na saída do Blast")
    print("Saiu da funcao:", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("###################\n") 



'''
Função que exibe os sítios conservativos
@param self classe de onde foi chamado
@param results lista com o resultado da pesquisa lido por SeqIO e objeto em formato de dicionário que retornou do alinhamento
@param seq_count número de sequências que aparecerão na tela
'''
def ShowSites(self, results, seq_count):
    print("###################")
    print("Entrou na função de exibir os sitios:", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)

    search_record = results[0]
    blast_record = results[1]

    count=0
    queryies=[]
    matches=[]
    subjects=[]
    lengths=[]
    q_starts=[]
    s_starts=[]
    leng = 0
    maxs = 0

    path = "results/sites_result.txt"


    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if count < seq_count:
                count+=1


                queryies.append(hsp.query)
                q_starts.append(hsp.query_start)
                matches.append(hsp.match)
                subjects.append(hsp.sbjct)
                s_starts.append(hsp.sbjct_start)
                lengths.append(alignment.length)
                leng += alignment.length 

                if alignment.length > maxs:
                    maxs = alignment.length


    self.sitesFrame.setMinimumSize(QtCore.QSize(20*maxs, 1200))


    try:
        with open(path,'w') as result_file:
            result_file.write("\n*** Sítios ***\n\n")

            '''''' 
            # Escreve a sequência referência (o for é para percorrer os nuceotídeos e dar um espaco depois de mostrar cada um)
            result_file.write("*** Comparações entre a Query e as Subjects *** \n")
            result_file.write("Sequência:\tNúmero de mutações:\tBases:\n")

            result_file.write("Query:\t\t\t\t\t")

            ''' Escreve cada caractere da sequência com um espaço '''
            for c in range(0,len(queryies[0])):
                result_file.write(str(queryies[0][c]) + " ")

            result_file.write("\n")
            ''''''


            '''
            AAATTTCCCGGG     QUERY
            ||||||||||||    ||||||||
            AAATTTCCCGGG    SUBJECT0
            ||||||||||||    ||||||||
               TTTCCCGGG    SUBJECT1
            ||||||||||||    ||||||||
                ...           ...
                            SUBJECTN

            ~ Para cada sequência 's' de subjects
            ~ Para cada base 'base' da query
                - Se o índice da query for maior/igual ao índice que inicia a subject 
                    - Se a base da query com o índice correspondente ao índice da subject for igual a base da subject
                        Coloca a base
                    - Se a base da query com o índice correspondente ao índice da subject for diferente da base da subject
                        Coloca a base com um indicador de mutação
                - Se o índice não for igual/maior, coloca espaço em branco
            '''


            for i in range(len(subjects)):    # Percorre os indices de subjects
                result_file.write("Subject "+ str(i) +":\t\t")

                '''
                comparison(seq1,seq2) >>> Retorna uma lista com as posições onde ocorreu mutação
                '''
                positions = comparison(queryies[0],subjects[i])
                # positions = comparison("AAAATTTTCCCCGGGG","AAACTTTTCCCAGGGG") -> 2

                result_file.write(str( len(positions) ) + "\t\t")

                for c in range(0,len(subjects[i])): # Percorre as bases de cada subject

                    if queryies[0][c] == subjects[i][c]:
                        result_file.write( str(subjects[i][c]) + " ")
                    else:                                       
                        result_file.write(str(subjects[i][c]) + "!") # Mutação

                result_file.write("\n")

            result_file.write("\n\n")



    except:
        traceback.print_exc()


    ### Colocando a arquivo na tela do programa ###
    try:
        with open(path,"r") as r:
            self.sitesResultLabel.setText(r.read())
    
    except:
        traceback.print_exc()
    ###############################################


    print("Saiu da funcao: ", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("###################\n")




'''
Função que exibe os gráficos
@param self classe de onde foi chamado
@param results lista com o resultado da pesquisa lido por SeqIO e objeto em formato de dicionário que retornou do alinhamento
@param seq_count número de sequências que aparecerão na tela
'''
def ShowGraph(self, results, seq_count):

    print("###################")
    print("Entrou na função de exibir gráfico: ", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)

    search_record = results[0]
    blast_record = results[1]

    count = 0
    bases = {'A','C','G','T'}   # Dicionário com a quantidade de bases 
    sequence = []              
    sequences = []

    self.result_ui.basesFrame.setMinimumSize(QtCore.QSize(400, seq_count*600 ))

    try:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if count< seq_count:
                    count+=1

                    # sequence.append(alignment.title)

                    s = Seq(hsp.sbjct)

                    ###### Contando bases ######
                    count_a = s.count('A')
                    count_c = s.count('C')
                    count_t = s.count('T')
                    count_g = s.count('G')

                    bases = { 
                            'A': count_a,
                            'C': count_c,
                            'T': count_t,
                            'G': count_g 
                            }
                    
                    sequences.append(bases)
                    ############################
            # sequences.append(sequence)

    except :
        traceback.print_exc()

    print(sequences)

    bas = [1,2,3,4]
    bases_count = []




    ####### Criando um gráfico por sequência ########
    gs = []
    i=0

    for sequence in sequences:

        for base in sequence:
            bases_count.append(sequence[base])


        g = self.graph()
        pen = pg.mkPen(color=(255, 0, 0))
        bg1 = pg.BarGraphItem(x=bas, height=bases_count, width=0.3, brush='r', pen=pen)
        self.result_ui.graphicsView.setTitle( "Subject "+ str(i) + "\nA - T - G - C")
        self.result_ui.verticalLayout.addWidget(self.result_ui.graphicsView)
        g.addItem(bg1)
        gs.append(g)   

        bases_count.clear()
        i+=1


    print("Saiu da funcao: ", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("###################\n")


'''
Função Distância de Hamming que calcula o total de mutações
@param d1 primeira sequência 
@param d2 segunda sequência
'''
def hamming(d1, d2):
    total = 0

    if d1 > d2:
        menor_dna = d2
    else: menor_dna = d1
    
    for i in range(len( menor_dna )):
        if d1[i] != d2[i]:
            total += 1
    print(total)
    return total


'''
Função que compara e retorna uma lista com as posições das mutações
@param genome1 primeira sequência
@param genome2 segunda sequência
'''
def comparison(genome1, genome2):
    positions = []
    if (len(genome1) == len(genome2)):
        for i in range(len(genome1)):
            if genome1[i] != genome2[i]:
                positions.append(i)
        return positions

    if (len(genome1) > len(genome2)):
        for i in range(len(genome2)):
            if genome1[i] != genome2[i]:
                positions.append(i)
        if len(positions) == 0:
            positions.append(len(genome1) - 1)
        return positions
        
    if (len(genome1) < len(genome2)):
        for i in range(len(genome1)):
            if (genome1[i] != genome2[i]):
                positions.append(i)
        if len(positions) == 0:
            positions.append(len(genome2) - 1)
        return positions


'''
~ Ainda não
Função que exibe a árvore filogenética
@param self classe resultScreen
@param 
'''
def ShowPhylo(self,result_ui):
    print("Mostrando árvore filogenética..")
    # result_ui.resultText.setPixmap(QtGui.QPixmap("C:/Users/fabri_000/Documents/_Pesquisas TCC/Bioinformática Python/gui-pyqt5/images/loading_dna.gif"))
    
    
    # Alinhamento soemnte com arquivo local por enquanto
    # tree = Phylo.read('alignment_result.xml', 'phyloxml')
    # print(tree)

    print("Mostrada a árvore filogenética.")