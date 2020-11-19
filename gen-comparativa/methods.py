import os, sys, traceback

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QApplication, QMainWindow, QMessageBox
from PyQt5.QtCore import *

from pyqtgraph import PlotWidget, plot
import pyqtgraph as pg 

import Bio
from Bio import Entrez, SeqIO, AlignIO, Phylo  
from Bio.Phylo import PhyloXMLIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq

from validate_email import validate_email

from datetime import datetime

import math
import numpy as np


# Checa se tem formato de email
def CheckEmail(email):  
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



# Janela de pesquisa de arquivo
def OpenDialogBox(self,method_ui):
    options = QFileDialog.Options()
    filename = QFileDialog.getOpenFileName(filter="Fasta files (*.fasta)", options=options)
    path = filename[0]
    method_ui.insertEdit.setText(path)

    # print(self.filename.read())

    if path != '':
        record = SeqIO.read(path,"fasta")
        return record.id
    else:
        return "None"



def PopMultipleMethods(self):
    pop = QMessageBox()
    pop.setWindowTitle("Erro")
    pop.setText("Mais de um método selecionado")
    pop.setIcon(QMessageBox.Critical)
    pop.exec_()



def PopSearchError(self):
    pop = QMessageBox()
    pop.setWindowTitle("Erro")
    pop.setText("Erro na pesquisa")
    pop.setIcon(QMessageBox.Critical)
    pop.exec_()



def IsValidSearch(self):
    # [Válido/Não válido , Web/Local]

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

    handle = Entrez.esearch(db="nucleotide", term = method_ui.searchEdit.text(), idtype="acc", retmax = 1) # Retorna um XML
    record = Entrez.read(handle) #lendo as infos geradas pela pesquisa # Converte XML para estrutura de dados python (dicionários)
    handle.close()

    return record
    


def LocalAlignment(id):

    print("Entrou no alinhamento de arquivo local: ",datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("ID de alinhamento: ", id)

    result_handle = NCBIWWW.qblast("blastn", "nt", id)
    path = "results/alignment_result.xml"

    try:
        with open(path,"w") as save_file: # Para automaticamente fechar o "save_file"
            save_file.write(result_handle.read())
            result_handle.close()
    
    except IOError:
        print("Erro na escrita do arquivo de alinhamento: " + path)


    #########
    result_handle = open("/results/alignment_result_test.xml","r")
    #########

    try:
        with open(path,"r") as result_handle:
            blast_record = NCBIXML.read(result_handle)

    except IOError:
        print("Erro na leitura do arquivo de alinhamento: " + path)

    #########
    print("blast record: ") 
    print(blast_record)
    #########


    print("Fim do alinhamento: " + datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second + "\n")
    return blast_record



def WebAlignment(result_search):
    # Realiza o alinhamento do método de pesquisa Web
    print("######################")
    print("Entrei no alinhamento da pesquisa web:", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("ID de alinhamento: ", result_search['IdList'][0])

    # result_hande = NCBIWWW.qblast("blastn", "nt",result_search['IdList'][0], alignments=2)

    # path = "results/alignment_result.xml"
    path_test = "results/alignment_result_test.xml"

    # try:
    #     with open(path,"w") as save_file:
    #         save_file.write(result_hande.read())
    #         # print(type(save_file))
    # except IOError:
    #     print("Erro na escrita do arquivo de alinhamento: " + path)


    # try:
    #     with open(path,"r") as result_handle:
    #         result_handle = open(path,"r")
    #         blast_record = NCBIXML.read(result_handle)

    # except IOError:
    #     # Fazer rotina de erro na leitura
    #     print("Erro na leitura do arquivo de alinhamento " + path)
    #     # return False

    ########### Test ###########

    try:
        with open(path_test,"r") as result_handle:
            # exec('blast_record = NCBIXML.read(result_handle)')
            blast_record = NCBIXML.read(result_handle)

    except IOError:
        # Fazer rotina de erro na leitura
        traceback.print_exc()

    #############################


    # Caso não seja salvo em um arquivo, retornar a variável
    print("Fim alinhamento: ", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("######################\n")

    return blast_record



def ShowAlignments(self, blast_record, seq_count):
    print("###################")
    print("Entrou na função de exibir o alinhamento:", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)

    count = 0
    maxs = 0
    queryies=[]
    matches=[]
    subjects=[]
    lengths = []


    ### Para escrever na tela, é preciso gravar em um arquivo ###
    path = "results/result.txt"

    ######### Teste #########
    # print(blast_record)
    # print(blast_record.alignments)

    # for key in blast_record.alignments:
    #     print(key)

    # print(len(blast_record.alignments))
    #########################

    #################
    # try:
    #     for alignment in blast_record.alignments:
    #         for hsp in alignment.hsps:
    #             print("Score: ", hsp.score)
    #             print("Bits: ", hsp.bits)
    #             print("E: ", hsp.expect)
    #             print("num_alignments: ", hsp.num_alignments)
    #             print("Identities: ", hsp.identities)
    #             print("Positives: ", hsp.positives)
    #             print("Gaps: ", hsp.gaps)
    #             print("Strand: ", hsp.strand)
    #             print("Frame: ",hsp.frame)
    #             print("Query start: ",hsp.query_start)
    #             print("Sbjct star: ", hsp.sbjct_start)
    #             print("\n\n")

    # except:
    #     traceback.print_exc()
    #     pass

    #################

    try:
        with open(path ,'w') as result_file:
            for alignment in blast_record.alignments:
                if count< seq_count:
                    for hsp in alignment.hsps:
                        count+=1


                        ########### Cabeçalho ###########
                        result_file.write("\n *** Alinhamento *** \n" +
                        "Sequência: " + str(alignment.title) + "\n" + 
                        "Comprimento: "+ str(alignment.length) + "\n"+
                        "Número de bases diferentes: " + str(hsp.gaps)+"\n\n")


                        ########### Sequências ###########
                        # Para escrever do mesmo jeito que no BLAST
                        # math.ceil() arredonda o número para cima
                        for b in range( 0, math.ceil( len(hsp.query) / 60 )):

                            # q_s = hsp.query_start

                            result_file.write(
                                str(hsp.query_start + ((60+1) * b ) ) + " " +  str(hsp.query[ b*60 : (b+1)*60 ] ) + " " + str((hsp.query_start + 60) * (b+1)) + "\n" + 
                                str(hsp.sbjct_start + ((60+1) * b ) ) + " " + str(hsp.sbjct[ (b*60) : (b+1)*60 ] ) + " " + str((hsp.sbjct_start + 60) * (b+1)) + "\n\n")

                        result_file.write("\n\n")

                        ##########
                        queryies.append(hsp.query)
                        matches.append(hsp.match)
                        subjects.append(hsp.sbjct)

                        if alignment.length > maxs:
                            maxs = alignment.length
                        ##########

    except IOError:
        print("Erro na escrita: " + path)


    self.alignmentScrollFrame.setMinimumSize(QtCore.QSize(1200,  1.5*maxs*(count/3) )) #count*15000)) 


    try:
        with open(path,"r") as r:
            self.alignmentResultLabel.setText(r.read())

    except IOError:
        print("Erro na leitura: " + path)


    # print("\n" + "Tem ",count,"sequências na saída do Blast")
    print("Saiu da funcao:", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("###################\n") 



def ShowSites(self, blast_record, seq_count):
    print("###################")
    print("Entrou na função de exibir os sitios:", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)

    count=0
    queryies=[]
    matches=[]
    subjects=[]
    lengths=[]
    leng = 0

    path = "results/sites_result.txt"


    #################
    # try:
    #     for alignment in blast_record.alignments:
    #         for hsp in alignment.hsps:
    #             print("Score: ", hsp.score)
    #             print("Bits: ", hsp.bits)
    #             print("E: ", hsp.expect)
    #             print("num_alignments: ", hsp.num_alignments)
    #             print("Identities: ", hsp.identities)
    #             print("Positives: ", hsp.positives)
    #             print("Gaps: ", hsp.gaps)
    #             print("Strand: ", hsp.strand)
    #             print("Frame: ",hsp.frame)
    #             print("Query start: ",hsp.query_start)
    #             print("Sbjct star: ", hsp.sbjct_start)
    #             print("\n\n")

    # except:
    #     traceback.print_exc()
    #     pass

    #################





    try:
        with open(path,'w') as result_file:
            result_file.write("\n *** Sítios *** \n\n")

            for alignment in blast_record.alignments:
                if count< seq_count:
                    for hsp in alignment.hsps:
                        count+=1

                        queryies.append(hsp.query)
                        matches.append(hsp.match)
                        subjects.append(hsp.sbjct)
                        lengths.append(alignment.length)
                        leng += alignment.length 


            max = 0;

            for i in range(0,len(lengths)):
                if lengths[i] > max:
                    max = lengths[i]

            # print(max)

            self.sitesFrame.setMinimumSize(QtCore.QSize(20*max, 1200))

            # Escreve a sequência referência (o for é para percorrer os nuceotídeos e dar um espaco depois de mostrar cada um)
            result_file.write("Query:\t\t")

            for c in range(0,len(queryies[0])):
                result_file.write(str(queryies[0][c]) + " ")

            result_file.write("\n")



            # Escreve os subjects (um para cada elemento do)
            for i in range(0,len(subjects)):    # Percorre os indices de subjects
                result_file.write("Subject "+ str(i) +":\t")
                for c in range(0,len(subjects[i])): # Percorre as bases de cada subject
                    if queryies[0][c] == subjects[i][c]:
                        result_file.write( str(subjects[i][c]) + " ")
                    else:   # Mutação
                            # Falta um jeito de marcar a letra sem que fique desalinhado
                            # Se trocar de cor, troca de todas outras letrar na label
                        result_file.write(str(subjects[i][c]) + "! ")

                result_file.write("\n")
    
    except:
        traceback.print_exc()


    try:
        with open(path,"r") as r:
            self.sitesResultLabel.setText(r.read())
    
    except:
        traceback.print_exc()


    print("Saiu da funcao: ", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("###################\n")


def ShowGraph(self, blast_record, seq_count):

    print("###################")
    print("Entrou na função de exibir gráfico: ", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)

    count = 0
    bases = {'A','T','G','C'}   # Dicionário com a quantidade de bases 
    sequences = []              # Vetor com os dicíonários -> [ números de bases da sequência 1 , números de bases da sequência 2 ...]

    self.result_ui.basesFrame.setMinimumSize(QtCore.QSize(400, seq_count*1000 ))

    try:
        for alignment in blast_record.alignments:
            if count< seq_count:
                for hsp in alignment.hsps:
                    count+=1

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

    except IOError:
        print("Erro no desenho dos gráficos")

    # print(sequences)

    bas = [1,2,3,4]
    bases_count = []


    for base in bases:
        bases_count.append(bases[base])


    ####### Criando um gráfico por sequência ########
    gs = []

    for sequence in sequences:
        g = self.graph()
        bg1 = pg.BarGraphItem(x=bas, height=bases_count, width=0.6, brush='r')
        self.result_ui.verticalLayout.addWidget(self.result_ui.graphicsView)
        g.addItem(bg1)
        gs.append(g)    


    print("Saiu da funcao: ", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    print("###################\n")



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


def ShowPhylo(self,result_ui):
    print("Mostrando árvore filogenética..")
    # result_ui.resultText.setPixmap(QtGui.QPixmap("C:/Users/fabri_000/Documents/_Pesquisas TCC/Bioinformática Python/gui-pyqt5/images/loading_dna.gif"))
    
    


    # Alinhamento soemnte com arquivo local por enquanto
    # tree = Phylo.read('alignment_result.xml', 'phyloxml')
    # print(tree)

    print("Mostrada a árvore filogenética.")