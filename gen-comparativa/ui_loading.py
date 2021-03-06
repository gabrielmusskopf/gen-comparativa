# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'loadingWindow_v01.ui'
#
# Created by: PyQt5 UI code generator 5.15.1
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_LoadingWindow(object):
    def setupUi(self, LoadingWindow):
        LoadingWindow.setObjectName("LoadingWindow")
        LoadingWindow.resize(559, 358)
        self.centralwidget = QtWidgets.QWidget(LoadingWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.back_frame = QtWidgets.QFrame(self.centralwidget)
        self.back_frame.setStyleSheet("QFrame{\n"
"    background-color: rgb(30, 30, 30);\n"
"}")
        self.back_frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.back_frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.back_frame.setObjectName("back_frame")
        self.gridLayout = QtWidgets.QGridLayout(self.back_frame)
        self.gridLayout.setObjectName("gridLayout")
        self.loadingFrame = QtWidgets.QFrame(self.back_frame)
        self.loadingFrame.setMaximumSize(QtCore.QSize(400, 250))
        self.loadingFrame.setStyleSheet("QFrame{\n"
"    background-color: rgb(40, 40, 40);\n"
"    border-radius: 5px;\n"
"}\n"
"")
        self.loadingFrame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.loadingFrame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.loadingFrame.setObjectName("loadingFrame")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.loadingFrame)
        self.verticalLayout_2.setContentsMargins(-1, -1, -1, 20)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.titleLabel = QtWidgets.QLabel(self.loadingFrame)
        self.titleLabel.setMaximumSize(QtCore.QSize(16777215, 50))
        font = QtGui.QFont()
        font.setFamily("Bahnschrift SemiBold")
        font.setPointSize(30)
        self.titleLabel.setFont(font)
        self.titleLabel.setStyleSheet("QLabel{\n"
"    \n"
"    color: rgb(254, 165, 223);\n"
"}")
        self.titleLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.titleLabel.setObjectName("titleLabel")
        self.verticalLayout_2.addWidget(self.titleLabel)
        self.descriptionLabel = QtWidgets.QLabel(self.loadingFrame)
        self.descriptionLabel.setMaximumSize(QtCore.QSize(16777215, 20))
        font = QtGui.QFont()
        font.setFamily("Segoe UI")
        font.setPointSize(12)
        self.descriptionLabel.setFont(font)
        self.descriptionLabel.setStyleSheet("color: rgb(98, 114, 164);")
        self.descriptionLabel.setTextFormat(QtCore.Qt.AutoText)
        self.descriptionLabel.setScaledContents(False)
        self.descriptionLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.descriptionLabel.setWordWrap(True)
        self.descriptionLabel.setIndent(-1)
        self.descriptionLabel.setObjectName("descriptionLabel")
        self.verticalLayout_2.addWidget(self.descriptionLabel)
        self.imageLabel = QtWidgets.QLabel(self.loadingFrame)
        font = QtGui.QFont()
        font.setPointSize(20)
        self.imageLabel.setFont(font)
        self.imageLabel.setStyleSheet("")
        self.imageLabel.setAlignment(QtCore.Qt.AlignCenter)
        self.imageLabel.setObjectName("imageLabel")
        self.verticalLayout_2.addWidget(self.imageLabel)
#         self.progressBar = QtWidgets.QProgressBar(self.loadingFrame)
#         self.progressBar.setStyleSheet("QProgressBar {\n"
# "    \n"
# "    background-color: rgb(98, 114, 164);\n"
# "    color: rgb(200, 200, 200);\n"
# "    border-style: none;\n"
# "    border-radius: 10px;\n"
# "    text-align: center;\n"
# "}\n"
# "QProgressBar::chunk{\n"
# "    border-radius: 10px;\n"
# "    background-color: qlineargradient(spread:pad, x1:0, y1:0.511364, x2:1, y2:0.523, stop:0 rgba(254, 121, 199, 255), stop:1 rgba(170, 85, 255, 255));\n"
# "}")
#         self.progressBar.setProperty("value", 24)
#         self.progressBar.setObjectName("progressBar")
#         self.verticalLayout_2.addWidget(self.progressBar)
        self.gridLayout.addWidget(self.loadingFrame, 0, 0, 1, 1)
        self.verticalLayout.addWidget(self.back_frame)
        LoadingWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(LoadingWindow)
        QtCore.QMetaObject.connectSlotsByName(LoadingWindow)

    def retranslateUi(self, LoadingWindow):
        _translate = QtCore.QCoreApplication.translate
        LoadingWindow.setWindowTitle(_translate("LoadingWindow", "MainWindow"))
        self.titleLabel.setText(_translate("LoadingWindow", "Carregando"))
        self.descriptionLabel.setText(_translate("LoadingWindow", "<strong>CARREGANDO </strong> ALINHAMENTOS"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    LoadingWindow = QtWidgets.QMainWindow()
    ui = Ui_LoadingWindow()
    ui.setupUi(LoadingWindow)
    LoadingWindow.show()
    sys.exit(app.exec_())
