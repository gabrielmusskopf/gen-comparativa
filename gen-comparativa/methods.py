import os

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *

import Bio
from Bio import Entrez, SeqIO, AlignIO, Phylo  
from Bio.Phylo import PhyloXMLIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalOmegaCommandline

from validate_email import validate_email

from datetime import datetime

import math


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

    print("Entrei no LocalAlignment: ",datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)

    # result_handle = NCBIWWW.qblast("blastn", "nt", id, "XML")

    # with open("/results/alignment_result.xml","w") as save_file: # Para automaticamente fechar o "save_file"
    #     save_file.write(result_handle.read())
    #     result_handle.close()


    #########
    result_handle = open("/results/alignment_result_test.xml","r")
    #########
    # result_handle = open("/results/alignment_result.xml","r")
    blast_record = NCBIXML.read(result_handle)

    #########
    print("blast record: ") 
    print(blast_record)
    #########


    print("Fim do LocalAlignment: ", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)
    return blast_record



def WebAlignment(result_search):
    # Realiza o alinhamento do método de pesquisa Web

    print("Entrei no WebAlignment:", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)

    # result_hande = NCBIWWW.qblast("blastn", "nt",result_search['IdList'][0], alignments=2)

    # save_file = open("/results/alignment_result.xml","w")
    # save_file.write(result_hande.read())
    # # print(type(save_file))
    # save_file.close()

    try:
        path = "results/alignment_result.xml"
        result_handle = open(path,"r")
        # result_handle = open("/results/alignment_result_test.xml","r")
        blast_record = NCBIXML.read(result_handle)
        ########
        print(blast_record)
        ########

    except IOError:
        # Fazer rotina de erro na leitura
        print("Problem reading: " + path)
        # return False

    


    # Caso não seja salvo em um arquivo, retornar a variável
    print("Fim do WebAlignment: ", datetime.now().hour, ":", datetime.now().minute, ":", datetime.now().second)

    return blast_record

'''
Falta trocar de tela quando o alinhamento é local
'''

def ShowAlignments(self, blast_record, methodScrn):

    count = 0
    toShow = int( methodScrn.method_ui.sequencesAskLineEdit_2.text() )
    maxs = 0
    queryies=[]
    matches=[]
    subjects=[]
    lengths = []

    # print( int(methodScrn.method_ui.sequencesAskLineEdit_2.text()) )

    path = "results/result.txt"

    try:
        with open(path ,'w') as result_file:
            for alignment in blast_record.alignments:
                if count< toShow:
                    for hsp in alignment.hsps:
                        count+=1

                        ########### Cabeçalho ###########
                        result_file.write("\n *** Alinhamento *** \n" +
                        "Sequência: " + str(alignment.title) + "\n" + 
                        "Comprimento: "+ str(alignment.length) + "\n"+
                        "Número de bases diferentes: "+str(hsp.gaps)+"\n")


                        ########### Sequências ###########
                        # Para escrever do mesmo jeito que no BLAST
                        # math.ceil() arredonda o número para cima
                        for b in range( 0, math.ceil( len(hsp.query) / 60 )):
                            result_file.write(
                                str((b*60)+1) + " " +  str(hsp.query[b*60:(b+1)*60]) + " " + str((b+1)*60) + "\n" + 
                                str((b*60)+1) + " " + str(hsp.sbjct[(b*60):(b+1)*60]) + " " + str((b+1)*60) + "\n\n")

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


    print("\n" + "Tem ",count,"sequências na saída do Blast")



def ShowSites(self,blast_record, methodScrn):

    count=0
    toShow = int(methodScrn.method_ui.sequencesAskLineEdit_2.text())
    queryies=[]
    matches=[]
    subjects=[]
    lengths=[]
    leng = 0

    path = "results/result.txt"

    try:
        with open(path,'w') as result_file:
            result_file.write("\n *** Sítios *** \n\n")

            for alignment in blast_record.alignments:
                if count< toShow:
                    for hsp in alignment.hsps:
                        count+=1
                        # self.result_ui.sitesFrame.setMinimumSize(QtCore.QSize(12*alignment.length, 1200))
                        # result_file.write( str(hsp.query) + "\t" +"\n")
                        # result_file.write( str(hsp.match) + "\t" +"\n")
                        queryies.append(hsp.query)
                        matches.append(hsp.match)
                        subjects.append(hsp.sbjct)
                        lengths.append(alignment.length)
                        leng += alignment.length 

            # queryies=['ATGCATGCATGCATGC','GATCGATCGATCGATC','ATGCTAGCTAGTCAT']
            # subjects=['ATGGATGCATGCTTGC','GATTGATCGCTCGATC','ATGATAGCTAGTCCT']
            # lengths=[len(queryies[0]),len(queryies[1]),len(queryies[2])]

            # print(lengths)

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
    
    except IOError:
        print("Erro na escrita: " + path)


    try:
        with open(path,"r") as r:
            self.sitesResultLabel.setText(r.read())
    
    except IOError:
        print("Erro na leitura: " + path)

    # Para árvore filogenética
    # wdir = os.getcwd()
    # print(wdir)

    # in_file = "sites_result.txt"
    # out_file = "alignment_general.fasta"

    # clustalomega_cline = ClustalOmegaCommandline(infile = in_file, outfile = out_file, verbose = True, auto = False)
    # print(clustalomega_cline)

    # path = "r'" + str(wdir)
    # input_comand = '"' + wdir + '/clustal-omega-1.2.2-win64/' + str(clustalomega_cline)[0:8] + '" ' + str(clustalomega_cline)
    # print(input_comand)
    # os.system(input_comand)



    # for i in range(0,len(queryies)):
    #     for j in range(lengths[0]):
    #         if queryies[0][j] != subjects[i][j]:
    #             print("Alteração na posição: " + str(i) + ", " + str(j) + "\n")

    
    # print(sequences[0][0])




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