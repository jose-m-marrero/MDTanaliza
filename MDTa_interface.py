#!/usr/bin/env python
#-*- coding: utf-8 -*-

import os
import sys
from Tkinter import *
from Tkinter import Tk, Frame, BOTH
import Tkconstants, tkFileDialog
from tkMessageBox import showerror, showinfo, askyesno
import numbers


"""
* Copyright (C) 2009-2017  Jose M. Marrero <josemarllin@gmail.com> and Ramon Ortiz <ramonortizramis@gmail.com>
* You may use, distribute and modify this code under the terms of the MIT License.
* The authors will not be held responsible for any damage or losses or for any implications 
* whatsoever resulting from downloading, copying, compiling, using or otherwise handling this 
* source code or the program compiled.
Code name: MDTanaliza GUI
Version: 2.1-2018-04-06
Execute command line:  ./MDTa_interface.py
Authors: Jose M. Marrero (1) R. Ortiz Ramis (2)
Affiliation 1: REPENSAR
Affiliation 1: IGEO-CSIC
"""

'''
x = StringVar() # Holds a string; default value ""
x = IntVar() # Holds an integer; default value 0
x = DoubleVar() # Holds a float; default value 0.0
x = BooleanVar() # Holds a boolean, returns 0 for False and 1 for True
To read the current value of such a variable, call the method get(). The value of such a variable can be changed with the set() method. 
'''

class Example(Frame):
	
	def __init__(self, parent):
		Frame.__init__(self, parent)   
		self.parent = parent
		self.ntipe = 0
		self.LAN = 2
		if self.LAN < 1 or self.LAN > 2:
			print "Atencion, el numero asignado al lenguage debe estar entre 1 o 2"
			sys.exit()
		self.initUI()   
	#CHAGE LANGUAGE
	def chalangen2sp(self):
		if (self.LAN == 2):
			self.LAN = 1
			self.initUI() 
			
	def chalangsp2en(self):		
		if (self.LAN == 1):
			self.LAN = 2
			self.initUI()	
			
	
	#UPDATE DISABLE-NORMAL----------------------------------------------	
	def update_chk(self):
		"""Change status according to weather some buttons are pressed or not."""
		#MODIFICA MDT
		if self.newz.get():
			self.fase1.config(state=NORMAL)
			self.fase2.config(state=NORMAL)
			if self.nfase.get() == 1:
				self.dirma.delete(0,END)
				self.dirma.config(state=NORMAL)
				self.nwz.delete(0,END)
				self.nwz.config(state=DISABLED)
				self.maskButton.config(state=NORMAL)
				self.newButton.config(state=DISABLED) 
				self.masktipe.set(0)
				self.maskradb.config(state=NORMAL)
				self.maskrada.config(state=NORMAL)
			if self.nfase.get() == 2:
				self.dirma.delete(0,END)
				self.dirma.config(state=DISABLED)
				self.nwz.delete(0,END)
				self.nwz.config(state=NORMAL)
				self.maskButton.config(state=DISABLED)
				self.newButton.config(state=NORMAL)
				self.masktipe.set(0)
				self.maskradb.config(state=DISABLED)
				self.maskrada.config(state=DISABLED)
		if self.newz.get() == 0:
			self.nfase.set(0)
			self.fase1.config(state=DISABLED)
			self.fase2.config(state=DISABLED)
			self.masktipe.set(0)
			self.maskradb.config(state=DISABLED)
			self.maskrada.config(state=DISABLED)
			self.dirma.delete(0,END)	
			self.dirma.config(state=DISABLED)
			self.maskButton.config(state=DISABLED)
			self.nwz.delete(0,END)
			self.nwz.config(state=DISABLED)
			self.newButton.config(state=DISABLED) 	
		#RECORTE
		if self.recor.get():
			self.xmin.config(state=NORMAL) 
			self.xmax.config(state=NORMAL) 
			self.ymin.config(state=NORMAL) 
			self.ymax.config(state=NORMAL)
			self.resx.config(state=NORMAL)
			self.resy.config(state=NORMAL)
			
		if (self.recor.get() == 0):
			self.xmin.delete(0,END)
			self.xmin.config(state=DISABLED) 
			self.xmax.delete(0,END)
			self.xmax.config(state=DISABLED) 
			self.ymin.delete(0,END)
			self.ymin.config(state=DISABLED) 
			self.ymax.delete(0,END)
			self.ymax.config(state=DISABLED)
			self.resx.delete(0,END)
			self.resx.config(state=DISABLED)
			self.resy.delete(0,END)
			self.resy.config(state=DISABLED)
		#PROCESADO
		#sink -----------------------------
		valsik = int(self.sikvar.get())
		if valsik == 0:
			self.siklab.config(bg="#d9d9d9", foreground="black", text=self.txtchk9)
		if valsik == 1:
			self.siklab.config(bg="yellow", foreground="blue", text=self.txtchk9a)
		if valsik == 2:
			self.siklab.config(bg="yellow", foreground="blue", text=self.txtchk9b)	
		#aspcet ----------------------------
		valasp = int(self.aspvar.get())
		if valasp == 0:
			self.asplab.config(bg="#d9d9d9", foreground="black", text=self.txtchk10)
		if valasp == 1:
			self.asplab.config(bg="yellow", foreground="blue", text=self.txtchk10a)
		if valasp == 2:
			self.asplab.config(bg="yellow", foreground="blue", text=self.txtchk10b)		
		#slope	-----------------------------
		valslo = int(self.slopvar.get())
		if valslo == 0:
			self.sloplab.config(bg="#d9d9d9", foreground="black", text=self.txtchk11)
		if valslo == 1:
			self.sloplab.config(bg="yellow", foreground="blue", text=self.txtchk11a)
		if valslo == 2:
			self.sloplab.config(bg="yellow", foreground="blue", text=self.txtchk11b)
		if valslo == 3:
			self.sloplab.config(bg="yellow", foreground="blue", text=self.txtchk11c)	
		if valslo == 4:
			self.sloplab.config(bg="yellow", foreground="blue", text=self.txtchk11d)
		if valslo == 5:
			self.sloplab.config(bg="yellow", foreground="blue", text=self.txtchk11e)
		if valslo == 6:
			self.sloplab.config(bg="yellow", foreground="blue", text=self.txtchk11f)
		if valslo == 7:
			self.sloplab.config(bg="yellow", foreground="blue", text=self.txtchk11g)
		if valslo == 8:
			self.sloplab.config(bg="yellow", foreground="blue", text=self.txtchk11h)	
		if valslo == 9:
			self.sloplab.config(bg="yellow", foreground="blue", text=self.txtchk11i)			
					
		#TRAYECTO
		#algor = int(self.npt.get()) 
		algor = int(self.tratip.get())
		if(algor == 0):
			self.dismax.delete(0,END)
			self.dismax.config(state=DISABLED)
			self.altmax.delete(0,END)
			self.altmax.config(state=DISABLED)
			self.rad.delete(0,END)
			self.rad.config(state=DISABLED)
			self.incre.delete(0,END)
			self.incre.config(state=DISABLED)
			self.force.set("0")                               #<-------------
			self.forceval.config(state=DISABLED)              #<-------------
			self.itera.delete(0,END)
			self.itera.config(state=DISABLED)
			self.npt.delete(0,END)
			self.npt.config(state=DISABLED)
			self.T.delete("0.0",'end')
			self.T.config(state=DISABLED)
			self.loadnpt.config(state=DISABLED)
			self.labalgtyp.config(text=self.txtlab18)
			self.mod.set("0")                               #<-------------
			self.modval.config(state=DISABLED)              #<-------------
			#self.mod.delete(0,END)
			#self.mod.config(state=DISABLED)
			self.huso.config(state=DISABLED)
			
			self.labhemis.config(bg="#d9d9d9", foreground="black", text=self.txtlab19) #<-------------
			self.hemis.set("0")	                                                       #<-------------
			self.hemisval.config(state=DISABLED)                                       #<-------------
			#self.hemis.config(state=DISABLED)                                         #<-------------
		#activado
		if(algor > 0):
			self.dismax.config(state=NORMAL)
			self.altmax.config(state=NORMAL)
			self.npt.config(state=NORMAL)
			self.T.config(state=NORMAL)
			self.loadnpt.config(state=NORMAL)
			#self.mod.config(state=NORMAL)
			self.modval.config(state="readonly")                        #<-------------
			self.huso.config(state=NORMAL)
			self.hemisval.config(state="readonly")
			
			typhem = int(self.hemis.get())                                                  #<-------------
			if(typhem == 0):                                                                #<-------------
				self.labhemis.config(bg="yellow", foreground="blue", text=self.txtlab19a)	#<-------------
			if(typhem == 1):                                                                #<-------------
				self.labhemis.config(bg="yellow", foreground="blue",text=self.txtlab19b)    #<-------------
			#self.hemis.config(state=NORMAL)                                                #<-------------
			
			if(algor == 1 or algor == 2):
				self.incre.config(state=NORMAL)
				self.forceval.config(state="readonly")            #<-------------
				self.rad.delete(0,END)
				self.rad.config(state=DISABLED)
				self.itera.delete(0,END)
				self.itera.config(state=DISABLED)
				if(algor == 1):
					self.labalgtyp.config(text=self.txtlab18a)
				if(algor == 2):
					self.labalgtyp.config(text=self.txtlab18b)		
			if(algor == 3):	
				self.incre.config(state=NORMAL)
				self.forceval.config(state="readonly")            #<-------------
				self.itera.config(state=NORMAL)	
				self.rad.delete(0,END)
				self.rad.config(state=DISABLED)	
				self.labalgtyp.config(text=self.txtlab18c)
			if(algor == 4):	
				self.rad.config(state=NORMAL)
				self.incre.delete(0,END)
				self.incre.config(state=DISABLED)
				self.force.set("0")                            #<-------------
				self.forceval.config(state=DISABLED)              #<-------------
				self.itera.delete(0,END)
				self.itera.config(state=DISABLED)
				self.labalgtyp.config(text=self.txtlab18d)
		
	#CHECK INTEGRITY----------------------------------------------------
	def checkseq(self):
		"""Secuencia de chequeo de informacion"""
		#check integrity sequence---------------------------------------
		self.totfall = 0
		#directorios
		dirdat = self.checkdirectories()
		#modifica
		if(self.newz.get() == 1):
			mod = self.checkmodifica()
		if(self.newz.get() == 0):
			msg2 = "Modifica no elegido -2-\n"		
			self.txtout.insert("0.0", msg2)
			print msg2
			mod = 1	
		#datos mdt	
		data = self.checkdatosdtm()
		#recorte
		if(self.recor.get() == 1):
			rec = self.checkrecorte()
		if(self.recor.get() == 0):
			msg2 = "Recorte no elegido -4-\n"		
			self.txtout.insert("0.0", msg2)
			print msg2
			rec = 1	
		#procesa	
		valsik = int(self.sikvar.get())
		valasp = int(self.aspvar.get())
		valslo = int(self.slopvar.get())
		if (valsik > 0 or valasp > 0 or valslo > 0):
			pros = self.checkprocesados()
		else:	
			msg2 = "Procesados no elegido -5-\n"		
			self.txtout.insert("0.0", msg2)
			print msg2
			pros = 1
		#trayec	
		if len(self.npt.get()) > 0:	
			val = int(self.npt.get())
			if(val > 0):
				tra = self.checktrayectorias()		
			if(val == 0):
				msg2 = "Trayectorias no elegidas -6-\n"		
				self.txtout.insert("0.0", msg2)
				print msg2	
				tra = 1	
		if len(self.npt.get()) == 0:
			msg2 = "Trayectorias no elegidas -6-\n"		
			self.txtout.insert("0.0", msg2)
			print msg2	
			tra = 1		
		glob = 	dirdat + mod + data + rec + pros + tra
		if glob == 6: #ok
			return 1
		if glob != 6: #errores
			return 0			
	
	def checkdirectories(self):
		self.demButton.configure(bg="#d9d9d9")
		self.outButton.configure(bg="#d9d9d9")
		fallos = 0
		#**********************************MDT
		self.valtx  = self.dirin.get()
		self.texlab = "MDT"
		if not self.dirin.get(): #no esta relleno
			self.is_filled()
			self.demButton.configure(bg="red")
			fallos +=1
		if self.dirin.get():	 #si esta relleno
			resul = self.is_filexist() #existe ruta
			if(resul == 0): #no existe
				self.demButton.configure(bg="red")
				fallos +=1
			if(resul == 1): #si existe
				if self.num < 2000:
					self.ntipe = self.is_binary() #es binario
					self.valnum = self.demtipe.get()
					if self.ntipe: #si es
						coin = self.binasc_ok() #coincide tipo
						if coin == 1: #no coincide
							self.labdtipe.configure(bg="red")
							self.demButton.configure(bg="red")
							fallos +=1	
					if not self.ntipe: #si no es
						coin = self.binasc_ok() #coincide tipo
						if coin == 2: #no coincide
							self.labdtipe.configure(bg="red")
							self.demButton.configure(bg="red")
							fallos +=1	
		#*************************
		self.valtx  = self.dirout.get()
		self.texlab = "Dirout"
		if not self.dirout.get(): #no esta relleno
			self.is_filled()
			self.outButton.configure(bg="red")
			fallos +=1
		if self.dirout.get():	 #si esta relleno
			resul = self.is_pathxist() #existe ruta
			if(resul == 0): #no existe
				self.outButton.configure(bg="red")
				fallos +=1				
		#****************
		if fallos > 0:
			msg2 = ("Rutas con %i fallos\n" % fallos)
			self.txtout.insert("0.0", msg2)
			print msg2
			self.totfall += fallos
			return 0
		if fallos == 0:
			msg2 = "Rutas correctas -1-\n"		
			self.txtout.insert("0.0", msg2)
			print msg2
			return 1
		
	def checkmodifica(self):
		self.maskButton.configure(bg="#d9d9d9")
		self.newButton.configure(bg="#d9d9d9")
		self.newbut.configure(bg="#d9d9d9")
		self.labmasktipe.configure(bg="#d9d9d9")
		#new xyz----------------------------------------------------
		fallos = 0
		if(self.newz.get() == 1):
			val = self.nfase.get()
			if val == 0:
				textmsg = "Atencion, no se ha elegido la fase"
				msg = showerror(title = "Error", message = textmsg)
				self.newbut.configure(bg="red")
				msg2 = "Error en archivo, falta el archivo XYZ\n"
				self.txtout.insert("0.0", msg2)
				print msg2
				fallos +=1	
			#****************************************mascara	
			if val == 1:	
				self.valtx  = self.dirma.get()
				self.texlab = "Mascara"
				if not self.dirma.get(): #no esta relleno
					self.is_filled()
					self.newbut.configure(bg="red")
					self.maskButton.configure(bg="red")
					fallos +=1
				if self.dirma.get():	 #si esta relleno
					resul = self.is_filexist() #existe ruta
					if(resul == 0): #no existe
						self.maskButton.configure(bg="red")
						fallos +=1
					if(resul == 1): #si existe
						self.ntipe = self.is_binary() #es binario
						if self.ntipe: #si es
							self.valnum = self.masktipe.get()
							coin = self.binasc_ok() #coincide tipo
							if coin == 1: #no coincide
								self.newbut.configure(bg="red")
								self.labmasktipe.configure(bg="red")
								fallos +=1	
						if not self.ntipe: #si no es
							self.valnum = self.masktipe.get()
							coin = self.binasc_ok() #coincide tipo
							if coin == 2: #no coincide
								self.newbut.configure(bg="red")
								self.labmasktipe.configure(bg="red")
								fallos +=1	
			#****************************************newxyz			
			if val == 2:
				self.valtx  = self.nwz.get()
				self.texlab = "NewXYZ"			
				if not self.nwz.get():
					self.is_filled()
					self.newbut.configure(bg="red")
					self.newButton.configure(bg="red")
					fallos +=1	
				if self.nwz.get():
					resul = self.is_filexist()
					if(resul == 0): #no existe
						self.newbut.configure(bg="red")
						self.newButton.configure(bg="red")
						fallos +=1						
				
			#****************	
			if fallos > 0:
				msg2 = ("Modifica con %i fallos\n" % fallos)
				self.txtout.insert("0.0", msg2)
				print msg2
				self.totfall += fallos
				return 0
			if fallos == 0:
				msg2 = "Modifica correcto -2-\n"		
				self.txtout.insert("0.0", msg2)
				print msg2
				return 1
		#**************	
		if(self.newz.get() == 0):
			msg2 = "Modifica no elegido -2-\n"		
			self.txtout.insert("0.0", msg2)
			print msg2
			return 1
						
	def checkdatosdtm(self):
		self.labdmax.configure(bg="#d9d9d9")
		self.labdmin.configure(bg="#d9d9d9")
		self.labdnull.configure(bg="#d9d9d9")
		#Valores DTM----------------------------------------------------
		fallos = 0
		#***********************************Maxval	
		totresp = 0
		self.valtx  = self.maxval.get()
		self.texlab = "Maxval"
		if len(self.maxval.get()) == 0:
			self.labdmax.configure(bg="red")
			self.is_filled()
			fallos +=1
		if len(self.maxval.get()) > 0:
			resp1 = self.is_punto()
			resp2 = self.is_decimal()
			totresp = resp1 + resp2
			if(totresp > 0):
				self.labdmax.configure(bg="red")
				fallos += totresp		
		#***********************************Minval	
		totresp = 0
		self.valtx  = self.minval.get()
		self.texlab = "Minimo"
		if len(self.minval.get()) == 0:
			self.labdmin.configure(bg="red")
			self.is_filled()
			fallos +=1
		if len(self.minval.get()) > 0:
			resp1 = self.is_punto()
			resp2 = self.is_decimal()
			totresp = resp1 + resp2
			if(totresp > 0):
				self.labdmin.configure(bg="red")
				fallos += totresp			
		#***********************************Nullval
		totresp = 0
		self.valtx  = self.nullval.get()
		self.texlab = "Nullval"
		if len(self.nullval.get()) == 0:
			self.labdnull.configure(bg="red")
			self.is_filled()
			fallos +=1
		if len(self.nullval.get()) > 0:
			resp1 = self.is_punto()
			resp2 = self.is_entero()
			if(resp1 == 1 or resp2 == 1):
				self.labdnull.configure(bg="red")
				if(resp1 == 1):
					fallos += 1
				if(resp2 == 1):
					fallos += 1	
			if(resp1 == 0 and resp2 == 0):
				self.valnum = int(self.nullval.get())
		#****************
		if fallos > 0:
			msg2 = ("MDT con %i fallos\n" % fallos)
			self.txtout.insert("0.0", msg2)
			print msg2
			self.totfall += fallos
			return 0
		if fallos == 0:
			msg2 = "Valores MDT correctos -3-\n"		
			self.txtout.insert("0.0", msg2)
			print msg2
			return 1
		
	def checkrecorte(self):	
		self.labdxmin.configure(bg="#d9d9d9")
		self.labdxmax.configure(bg="#d9d9d9")
		self.labdymin.configure(bg="#d9d9d9")
		self.labdymax.configure(bg="#d9d9d9")
		self.labdrex.configure(bg="#d9d9d9")
		self.labdrey.configure(bg="#d9d9d9")
		#recorte--------------------------------------------------------
		fallos = 0	
		if(self.recor.get() == 1):
			#************************************Xmin
			totresp = 0
			self.valtx  = self.xmin.get()
			self.texlab = "Xmin"
			if len(self.xmin.get()) == 0:
				self.labdxmin.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.xmin.get()) > 0:
				resp3 = 0
				resp1 = self.is_punto()
				resp2 = self.is_decimal()
				if((resp1 + resp2) == 0):
					self.valnum = float(self.xmin.get())
					resp3       = self.is_positivo()
				totresp = resp1 + resp2 + resp3
				if(totresp > 0):
					self.labdxmin.configure(bg="red")
					fallos += totresp	
			#************************************Xmax
			totresp = 0
			self.valtx  = self.xmax.get()
			self.texlab = "Xmax"
			if len(self.xmax.get()) == 0:
				self.labdxmax.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.xmax.get()) > 0:
				resp3 = 0
				resp1 = self.is_punto()
				resp2 = self.is_decimal()
				if((resp1 + resp2) == 0):
					self.valnum = float(self.xmax.get())
					resp3       = self.is_positivo()
				totresp = resp1 + resp2 + resp3
				if(totresp > 0):
					self.labdxmax.configure(bg="red")
					fallos += totresp			
			#************************************Ymin
			totresp = 0
			self.valtx  = self.ymin.get()
			self.texlab = "Ymin"
			if len(self.ymin.get()) == 0:
				self.labdymin.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.ymin.get()) > 0:
				resp3 = 0
				resp1 = self.is_punto()
				resp2 = self.is_decimal()
				if((resp1 + resp2) == 0):
					self.valnum = float(self.ymin.get())
					resp3       = self.is_positivo()
				totresp = resp1 + resp2 + resp3
				if(totresp > 0):
					self.labdymin.configure(bg="red")
					fallos += totresp		
			#*************************************Ymax
			totresp = 0
			self.valtx  = self.ymax.get()
			self.texlab = "Ymax"
			if len(self.ymax.get()) == 0:
				self.labdymax.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.ymax.get()) > 0:
				resp3 = 0
				resp1 = self.is_punto()
				resp2 = self.is_decimal()
				if((resp1 + resp2) == 0):
					self.valnum = float(self.ymax.get())
					resp3       = self.is_positivo()
				totresp = resp1 + resp2 + resp3
				if(totresp > 0):
					self.labdymax.configure(bg="red")
					fallos += totresp	
			#***************************************resx
			totresp = 0
			self.valtx  = self.resx.get()
			self.texlab = "Resx"
			if len(self.resx.get()) == 0:
				self.labdrex.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.resx.get()) > 0:
				resp3 = 0
				resp1 = self.is_punto()
				resp2 = self.is_decimal()
				if((resp1 + resp2) == 0):
					self.valnum = float(self.resx.get())
					resp3       = self.is_positivo()
				totresp = resp1 + resp2 + resp3
				if(totresp > 0):
					self.labdrex.configure(bg="red")
					fallos += totresp		
			#****************************************resy
			totresp = 0
			self.valtx  = self.resy.get()
			self.texlab = "Resy"
			if len(self.resy.get()) == 0:
				self.labdrey.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.resy.get()) > 0:
				resp3 = 0
				resp1 = self.is_punto()
				resp2 = self.is_decimal()
				if((resp1 + resp2) == 0):
					self.valnum = float(self.resy.get())
					resp3       = self.is_positivo()
				totresp = resp1 + resp2 + resp3
				if(totresp > 0):
					self.labdrey.configure(bg="red")
					fallos += totresp
			#****************	
			if fallos > 0:
				msg2 = ("Recorte con %i fallos\n" % fallos)
				self.txtout.insert("0.0", msg2)
				print msg2
				self.totfall += fallos
				return 0
			if fallos == 0:
				msg2 = "Recorte correcto -4-\n"		
				self.txtout.insert("0.0", msg2)
				print msg2
				return 1
		'''
		if(self.recor.get() == 0):
			msg2 = "Recorte no elegido -4-\n"		
			self.txtout.insert("0.0", msg2)
			print msg2	
			return 1	
		'''		
			
	def checkprocesados(self):
		#PROCESADO DTM----------------------------------------------
		fallos = 0
		used   = 0
		#****************************************sink
		self.valtx  = self.sikvar.get()
		self.texlab = self.txtchk9
		valsik = int(self.sikvar.get())	
		if valsik > 0:
			used = 1
			if valsik < 0 or valsik >2:
				textmsg = "Atencion, el valor debe estar entre 1-2"
				self.siklab.configure(bg="red")
				msg = showerror(title = "Error", message = textmsg)
				msg2 = ("Error, valor fuera de rango en %s\n" % self.txtchk9)
				self.txtout.insert("0.0", msg2)
				print msg2
				fallos +=1	
		#****************************************aspect
		self.valtx  = self.aspvar.get()
		self.texlab = self.txtchk10
		valasp = int(self.aspvar.get())
		if valasp > 0:
			used = 1
			if valasp < 0 or valasp >2:
				textmsg = "Atencion, el valor debe estar entre 1-2"
				self.asplab.configure(bg="red")
				msg = showerror(title = "Error", message = textmsg)
				msg2 = ("Error, valor fuera de rango en %s\n" % self.txtchk10)
				self.txtout.insert("0.0", msg2)
				print msg2
				fallos +=1	
		#****************************************slope
		self.valtx  = self.slopval.get()
		self.texlab = self.txtchk11	
		valslo = int(self.slopvar.get())
		if valslo > 0:
			used = 1
			if valslo < 0 or valslo >9:
				textmsg = "Atencion, el valor debe estar entre 1-9"
				self.sloplab.configure(bg="red")
				msg = showerror(title = "Error", message = textmsg)
				msg2 = ("Error, valor fuera de rango en slope\n" % self.txtchk11)
				self.txtout.insert("0.0", msg2)
				print msg2
				fallos +=1	
		#****************
		if used == 1:	
			if fallos > 0:
				msg2 = ("Procesados con %i fallos\n" % fallos)
				self.txtout.insert("0.0", msg2)
				self.totfall += fallos
				print msg2
				return 0
			if fallos == 0:
				msg2 = "Procesados correctos -5-\n"		
				self.txtout.insert("0.0", msg2)
				print msg2
				return 1
		if used == 0:		
			if (self.aspec.get() == 0):
				msg2 = "Procesados no elegido -5-\n"		
				self.txtout.insert("0.0", msg2)
				print msg2	
				return 1		
			
	def checktrayectorias(self):
		"""check integrity in the GUI"""		
		wkdir = os.getcwd()
		self.labdminlabdmax.configure(bg="#d9d9d9")
		self.labdminlabalt.configure(bg="#d9d9d9")
		self.labdminlabrad.configure(bg="#d9d9d9")
		self.labdminlabinc.configure(bg="#d9d9d9")
		self.labdminlabfor.configure(bg="#d9d9d9")       #<-------------
		self.labdminlabnpt.configure(bg="#d9d9d9")
		self.labdminlabite.configure(bg="#d9d9d9")
		self.T.configure(bg="white")
		
		#TRAYECTORIAS DTM-------------------------------------------
		fallos  = 0
		elegido = 0		
		algor = int(self.tratip.get()) #tipo algoritmo
		#Para cualquier situacion donde se active la opcion trayect
		if(algor > 0):	
			#************************************DistMax
			elegido = 1
			self.valtx  = self.dismax.get()
			self.texlab = "Dismax"
			if len(self.dismax.get()) == 0:
				self.labdminlabdmax.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.dismax.get()) > 0:
				resp3 = 0
				resp1 = self.is_punto()
				resp2 = self.is_decimal()
				if((resp1 + resp2) == 0):
					self.valnum = float(self.dismax.get())
					resp3       = self.is_positivo()
				totresp = resp1 + resp2 + resp3
				if(totresp > 0):
					self.labdminlabdmax.configure(bg="red")
					fallos += totresp	
			#************************************AlstMax
			self.valtx  = self.altmax.get()
			self.texlab = "Altmax"
			if len(self.altmax.get()) == 0:
				self.labdminlabalt.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.altmax.get()) > 0:
				resp3 = 0
				resp1 = self.is_punto()
				resp2 = self.is_decimal()
				if((resp1 + resp2) == 0):
					self.valnum = float(self.altmax.get())
					resp3       = self.is_positivo()
				totresp = resp1 + resp2 + resp3
				if(totresp > 0):
					self.labdminlabalt.configure(bg="red")
					fallos += totresp
			#************************************npt
			if len(self.npt.get()) == 0:
				textmsg = "Atencion, faltan los puntos para evaluar trayectorias"
				self.labdminlabnpt.configure(bg="red")
				msg = showerror(title = "Error", message = textmsg)
				msg2 = "Error en trayectorias, falta el algoritmo\n"
				self.txtout.insert("0.0", msg2)
				print msg2
				fallos += 1
			#***********************************XYZpt	
			if len(self.npt.get()) > 0:	
				puntos = int(self.npt.get())      #puntos
				content = self.T.get("1.0", END)
				totlineas = int(self.T.index('end-1c').split('.')[0]) 
				if len(content) <= 1:
					textmsg = "Atencion, faltan los puntos con coordenadas\n"
					self.T.configure(bg="red")
					msg = showerror(title = "Error", message = textmsg)
					msg2 = "Error en trayectorias, falta los puntos con coordenadas\n"
					self.txtout.insert("0.0", msg2)
					print msg2
					fallos += 1
				if totlineas < puntos:
					textmsg = "Atencion, el numero de coordenadas es menor que el indicado\n"
					self.labdminlabnpt.configure(bg="red")
					self.T.configure(bg="red")
					msg = showerror(title = "Error", message = textmsg)
					msg2 = "Error en trayectorias, faltan puntos con coordenadas\n"
					self.txtout.insert("0.0", msg2)
					print msg2
					fallos += 1		
			#***********************************MODE
			self.valtx  = self.mod.get()
			self.texlab = "Mode"
			if len(self.mod.get()) == 0:
				self.labdminlabmod.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.mod.get()) > 0:
				resp1 = self.is_punto()
				resp2 = self.is_entero()
				if(resp1 == 1 or resp2 == 1):
					self.labdminlabmod.configure(bg="red")
					if(resp1 == 1):
						fallos += 1
					if(resp2 == 1):
						fallos += 1
				var = int(self.mod.get())
				if var < 0 or var > 1:
					textmsg = "Atencion, MODE debe estar entre 0 y 1\n"
					self.labdminlabmod.configure(bg="red")
					msg = showerror(title = "Error", message = textmsg)
					msg2 = "Error en mode\n"
					self.txtout.insert("0.0", msg2)
					print msg2
					fallos += 1		
			#***********************************HUSO
			self.valtx  = self.huso.get()
			self.texlab = "Huso"
			if len(self.huso.get()) == 0:
				self.labhuso.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.huso.get()) > 0:
				resp1 = self.is_punto()
				resp2 = self.is_entero()
				if(resp1 == 1 or resp2 == 1):
					self.labhuso.configure(bg="red")
					if(resp1 == 1):
						fallos += 1
					if(resp2 == 1):
						fallos += 1
				var = int(self.huso.get())
				if var < 0 or var > 60:
					textmsg = "Atencion, el uso debe estar entre 0 y 60\n"
					self.labhuso.configure(bg="red")
					msg = showerror(title = "Error", message = textmsg)
					msg2 = "Error en huso\n"
					self.txtout.insert("0.0", msg2)
					print msg2
					fallos += 1
			#***********************************hemis
			self.valtx  = self.hemis.get()
			self.texlab = "Hemisferio"
			if len(self.hemis.get()) == 0:
				self.labhemis.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.hemis.get()) > 0:
				resp1 = self.is_punto()
				resp2 = self.is_entero()
				if(resp1 == 1 or resp2 == 1):
					self.labhemis.configure(bg="red")
					if(resp1 == 1):
						fallos += 1
					if(resp2 == 1):
						fallos += 1
				var = int(self.hemis.get())
				if var < 0 or var > 1:
					textmsg = "Atencion, el hemisferio es 0 norte y 1 sur\n"
					self.labhemis.configure(bg="red")
					msg = showerror(title = "Error", message = textmsg)
					msg2 = "Error en huso\n"
					self.txtout.insert("0.0", msg2)
					print msg2
					fallos += 1					
		#************************************RadBus
		if(algor == 4):
			self.valtx  = self.rad.get()
			self.texlab = "RadBus-Rest"
			if len(self.rad.get()) == 0:
				self.labdminlabrad.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.rad.get()) > 0:
				resp1 = self.is_punto()
				resp2 = self.is_decimal()
				if(resp1 == 1 or resp2 == 1):
					self.labdminlabrad.configure(bg="red")
					if(resp1 == 1):
						fallos += 1
					if(resp2 == 1):
						fallos += 1	
		#************************************IncAlt
		if(algor == 1 or algor == 2 or algor == 3):
			self.valtx  = self.incre.get()
			self.texlab = "IncAlt"
			if len(self.incre.get()) == 0:
				self.labdminlabinc.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.incre.get()) > 0:
				resp3 = 0
				resp1 = self.is_punto()
				resp2 = self.is_decimal()
				if((resp1 + resp2) == 0):
					self.valnum = float(self.incre.get())
					resp3       = self.is_positivo()
				totresp = resp1 + resp2 + resp3
				if(totresp > 0):
					self.labdminlabinc.configure(bg="red")
					fallos += totresp
					
		#***********************************FORCE
		if(algor == 1 or algor == 2 or algor == 3):
			self.valtx  = self.force.get()
			self.texlab = "Force"
			if len(self.force.get()) == 0:
				self.labdminlabfor.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.force.get()) > 0:
				resp1 = self.is_punto()
				resp2 = self.is_entero()
				if(resp1 == 1 or resp2 == 1):
					self.labdminlabfor.configure(bg="red")
					if(resp1 == 1):
						fallos += 1
					if(resp2 == 1):
						fallos += 1
				var = int(self.force.get())
				print var
				if var < 0 or var > 1:
					textmsg = "Atencion, FORCE debe estar entre 0 y 1\n"
					self.labdminlabfor.configure(bg="red")
					msg = showerror(title = "Error", message = textmsg)
					msg2 = "Error en mode\n"
					self.txtout.insert("0.0", msg2)
					print msg2
					fallos += 1					
		#************************************Itera
		if(algor == 3):	
			self.valtx  = self.itera.get()
			self.texlab = "Itera"
			if len(self.itera.get()) == 0:
				self.labdminlabite.configure(bg="red")
				self.is_filled()
				fallos +=1
			if len(self.rad.get()) > 0:
				resp1 = self.is_punto()
				resp2 = self.is_entero()
				if(resp1 == 1 or resp2 == 1):
					self.labdminlabite.configure(bg="red")
					if(resp1 == 1):
						fallos += 1
					if(resp2 == 1):
						fallos += 1				
		#************************************
		if(elegido == 1):								
			if fallos > 0:
				msg2 = ("Encontrados %i fallos\n" % fallos)
				self.txtout.insert("0.0", msg2)
				print msg2
				self.totfall += fallos
				return 0
			if fallos == 0:
				msg2 = "Trayectorias correctas -6-\n"		
				self.txtout.insert("0.0", msg2)
				print msg2
				return 1	
		if(elegido == 0):
			msg2 = "Trayectorias no elegidas -6-\n"		
			self.txtout.insert("0.0", msg2)
			print msg2	
			return 1	
	
	#CHECK BINARY - TYPE NUMBER-----------------------------------------
	def is_binary(self):
		filename = self.valtx
		fin = open(filename, 'rb')
		try:
			CHUNKSIZE = 1024
			while 1:
				chunk = fin.read(CHUNKSIZE)
				if '\0' in chunk: # found null byte
					return True #es binario
				if len(chunk) < CHUNKSIZE:
					break # done
		finally:
			fin.close()
		return False #es ascii
	
	def binasc_ok(self):
		
		formerror = 0
		if self.ntipe: #es bin
			if (self.valnum != 1 ):
				formerror = 1
		if not self.ntipe: #es asc
			if (self.valnum != 2):
				formerror = 1
		
		if	formerror == 1:	
			textmsg = "Atencion, el tipo (bin-ascii) no coincide con el archivo mascara"
			msg = showerror(title = "Error", message = textmsg)
			msg2 = "Error en tipo de archivo\n"
			self.txtout.insert("0.0", msg2)
			print msg2
			return 1
		if	formerror == 0:	
			return 0
				
	def is_filled(self):
		textmsg = ("Atencion, falta el valor %s" % self.texlab)
		msg = showerror(title = "Error", message = textmsg)
		msg2 = "Error en recorte, falta valor ymax\n"
		self.txtout.insert("0.0", msg2)
		print msg2
		return 1
		
	def is_punto(self):
		if ',' in self.valtx:
			textmsg = "Atencion, debe utilizar un punto para el separador decimal"
			msg = showerror(title = "Error", message = textmsg)
			msg2 = "Error en valor, punto para decimal\n"
			self.txtout.insert("0.0", msg2)
			print msg2
			return 1
		if not ',' in self.valtx:
			return 0	
	
	def is_decimal(self):
		if not '.' in self.valtx:	
			if self.ndectype == 0:
				textmsg = ("Atencion, el valor %s debe ser decimal" % self.texlab)
				msg = showerror(title = "Error", message = textmsg)
				msg2 = "Error en valor, debe ser decimal\n"
				self.txtout.insert("0.0", msg2)
				print msg2
				return 1
			if self.ndectype == 1:
				msg2 = "Archivo con cabecera\n"
				self.txtout.insert("0.0", msg2)
				self.ndectype = 0
				return 1	
		if '.' in self.valtx:
			self.ndectype = 0	
			return 0
	
	def is_entero(self):
		if '.' in self.valtx:
			if self.ndectype == 0:	
				textmsg = ("Atencion, el valor %s debe ser entero" % self.texlab)
				msg = showerror(title = "Error", message = textmsg)
				msg2 = "Error en valor, debe ser entero\n"
				self.txtout.insert("0.0", msg2)
				print msg2
				return 1
			if self.ndectype == 1:
				self.ndectype = 0
				return 1	
		if not '.' in self.valtx:
			self.ndectype = 0	
			return 0		
			
	def is_positivo(self):	
		if(self.valnum < 0):
			textmsg = ("Atencion, %s debe ser positivo" % self.texlab)
			msg = showerror(title = "Error", message = textmsg)
			msg2 = "Error en valor, debe ser postivio\n"
			self.txtout.insert("0.0", msg2)
			print msg2
			return 1
		if(self.valnum >= 0):
			return 0
	
	def is_negativo(self):	
		if(self.valnum >= 0):
			textmsg = ("Atencion, %s debe ser positivo" % self.texlab)
			msg = showerror(title = "Error", message = textmsg)
			msg2 = "Error en valor, debe ser postivio\n"
			self.txtout.insert("0.0", msg2)
			print msg2
			return 1
		if(self.valnum < 0):
			return 0				
	
	def is_filexist(self):
		filename = self.valtx
		#si es archivo
		if not os.path.isfile(filename):
			textmsg = ("Atencion, el archivo o la ruta a %s no existe" % self.texlab)
			msg = showerror(title = "Error", message = textmsg)
			msg2 = "Error en ruta al archivo\n"
			self.txtout.insert("0.0", msg2)
			print msg2
			return 0
		if os.path.isfile(filename):
			file_info = os.stat(filename)
			self.num = file_info.st_size
			#self.convert_bytes()
			return 1
	
	def is_pathxist(self):
		filename = self.valtx
		self.color = 0
		#si es archivo
		if not os.path.exists(filename):
			textmsg = (u"Atencion, la ruta a %s no existe, ¿desea crearla?" % self.texlab)
			msg = askyesno(title = "Error", message = textmsg)
			if msg == True:
				try:
					os.makedirs(filename)
					return 1
				except OSError as exception:
					if exception.errno != errno.EEXIST:
						return 0
			if msg == False:
				textmsg = (u"La ruta a %s no a sido creada" % self.texlab)
				msg = showerror(title = "Error", message = textmsg)
				msg2 = "La ruta no a sido creada\n"
				self.txtout.insert("0.0", msg2)
				print msg2
				return 0
		if os.path.exists(filename):
			return 1	
	
	def convert_bytes(self):
		"""
		this function will convert bytes to MB.... GB... etc
		"""
		for x in ['bytes', 'KB', 'MB', 'GB', 'TB']:
			if self.num < 1024.0:
				return "%3.1f %s" % (self.num, x)
			self.num /= 1024.0
					
	#CLEAN GUI----------------------------------------------------------	
	def cleangui(self):
		
		self.demButton.configure(bg="#d9d9d9")
		self.outButton.configure(bg="#d9d9d9")
		self.maskButton.configure(bg="#d9d9d9")
		self.newButton.configure(bg="#d9d9d9")
		self.newbut.configure(bg="#d9d9d9")
		self.labdmax.configure(bg="#d9d9d9")
		self.labdmin.configure(bg="#d9d9d9")
		self.labdnull.configure(bg="#d9d9d9")
		self.labdxmin.configure(bg="#d9d9d9")
		self.labdxmax.configure(bg="#d9d9d9")
		self.labdymin.configure(bg="#d9d9d9")
		self.labdymax.configure(bg="#d9d9d9")
		self.labdrex.configure(bg="#d9d9d9")
		self.labdrey.configure(bg="#d9d9d9")
		#---
		self.siklab.config(bg="#d9d9d9", foreground="black", text=self.txtchk9)
		self.asplab.config(bg="#d9d9d9", foreground="black", text=self.txtchk10)
		self.sloplab.config(bg="#d9d9d9", foreground="black", text=self.txtchk11)
		#---
		#self.labdminlabtray.configure(bg="#d9d9d9")
		self.labdminlabdmax.configure(bg="#d9d9d9")
		self.labdminlabalt.configure(bg="#d9d9d9")
		self.labdminlabrad.configure(bg="#d9d9d9")
		self.labdminlabinc.configure(bg="#d9d9d9")
		self.labdminlabmod.configure(bg="#d9d9d9")
		self.labdminlabnpt.configure(bg="#d9d9d9")
		self.T.configure(bg="white")
		self.labhuso.configure(bg="#d9d9d9")
		self.labhemis.configure(bg="#d9d9d9")
		
		self.dirin.delete(0,END)
		self.dirout.delete(0,END)
		
		self.newz.set(0)
		self.nfase.set(0)
		self.dirma.delete(0,END)
		self.dirma.config(state=DISABLED)	
		self.nwz.delete(0,END)
		self.nwz.config(state=DISABLED)
		self.newButton.config(state=DISABLED) 
		self.maskButton.config(state=DISABLED)	
		
		self.demtipe.set(0)
		#self.demradb.config(state=DISABLED)	
		#self.demrada.config(state=DISABLED)	
		self.maxval.delete(0,END)
		self.minval.delete(0,END)
		self.nullval.delete(0,END)
		
		self.recor.set(0)
		self.xmin.delete(0,END)
		self.xmin.config(state=DISABLED) 
		self.xmax.delete(0,END)
		self.xmax.config(state=DISABLED) 
		self.ymin.delete(0,END)
		self.ymin.config(state=DISABLED) 
		self.ymax.delete(0,END)
		self.ymax.config(state=DISABLED)
		self.resx.delete(0,END)
		self.resx.config(state=DISABLED)
		self.resy.delete(0,END)
		self.resy.config(state=DISABLED)
		
		self.sikvar.set('0')
		self.aspvar.set('0')
		self.slopvar.set('0')
		
		self.tratip.set(0)
		self.dismax.delete(0,END)
		self.altmax.delete(0,END)
		self.rad.delete(0,END)
		self.force.set("0") 
		self.incre.delete(0,END)
		self.mod.set('0')
		self.npt.delete(0,END)
		self.npt.insert( INSERT, "0" )
		self.T.delete("0.0",'end')
		self.huso.delete(0,END)
		self.huso.config(state=DISABLED)
		self.hemisval.delete(0,END)
		self.hemisval.config(state=DISABLED)
		
	
	#OPEN FILE----------------------------------------------------------
	def dtmfile(self):
		"""Returns a selected directoryname."""
		self.demButton.configure(bg="#d9d9d9")
		wkdir = os.getcwd()
		#txtdir = tkFileDialog.askdirectory(parent=self,title='Directorio entrada DTM')
		#comdir = txtdir + "/"
		filename = tkFileDialog.askopenfilename(parent=self,title='Nombre MDT')	
		if filename:
			self.txtout.insert("0.0", filename)
			#print os.path.commonprefix([wkdir, filename])
			finalfile = '/' + os.path.relpath(filename, wkdir)
			print finalfile
			print filename
			self.valtx = filename 
			self.dirin.delete(0,END)
			self.dirin.insert(INSERT, filename)
			self.ntipe = 1
			tipe = self.is_binary()
			if tipe:
				self.demtipe.set(1)
			else:
				self.demtipe.set(2)	
		if not filename:
			print "Accion cancelada\n"		
			
	def direout(self):
		"""Returns a selected directoryname."""
		self.outButton.configure(bg="#d9d9d9")
		txtdir = tkFileDialog.askdirectory(parent=self,title='Directorio de Salida')
		if txtdir:
			comdir = txtdir + "/"
			self.txtout.insert("0.0", comdir)
			self.dirout.delete(0,END)
			self.dirout.insert(INSERT, comdir)	
		if not txtdir:
			print "Accion cancelada\n"		
				
	def maskfile(self):
		"""Returns a selected directoryname."""
		self.maskButton.configure(bg="#d9d9d9")
		self.labmasktipe.configure(bg="#d9d9d9")
		#txtdir = tkFileDialog.askdirectory(parent=self,title='Directorio entrada DTM')
		#comdir = txtdir + "/"
		filename = tkFileDialog.askopenfilename(parent=self,title='Nombre Mascara')	
		if filename:
			self.txtout.insert("0.0", filename)
			print filename
			self.dirma.delete(0,END)
			self.dirma.insert(INSERT, filename)	
			self.ntipe = 2
			self.valtx = filename
			tipe = self.is_binary()
			if tipe:
				self.masktipe.set(1)
			else:
				self.masktipe.set(2)
		if not filename:
			print "Accion cancelada\n"			
		
	def newopenfile(self):	
		"""Returns a new xyz file name to modify the DEM."""
		self.newButton.configure(bg="#d9d9d9")
		self.newbut.configure(bg="#d9d9d9")
		filename = tkFileDialog.askopenfilename(parent=self,title='Archivo xyz nuevo')	
		if filename:
			print filename
			#name = os.path.basename(filename)
			self.txtout.insert("0.0", filename) 
			self.nwz.delete(0,END)
			self.nwz.insert(INSERT, filename)
		if not filename:
			print "Accion cancelada\n"		
	
	def openhelp(self):
		ostype = sys.platform
		if ostype == 'linux2':
			ejecutac = 'evince MDTa_help.pdf'
		if ostype == 'darwin':
			"open -a Preview.app MDTa_help.pdf"
		if ostype == "win32":
			print "windows"
		print ostype	
		os.system(ejecutac)	
		return 1
	
	#SAVE LOAD FILE-----------------------------------------------------
	def savefile(self):
		"""Save file, create a new file."""
		#check integrity--------------------------------------------
		valchk = self.checkseq()
		print(valchk)
		if valchk == 1:
			#escribe un fichero de salida
			self.cfgfile = tkFileDialog.asksaveasfilename(parent=self,filetypes=[('Config file','*.cfg')] ,title="Save configfile as...")		
			if self.cfgfile:
				self.txtfile =open(self.cfgfile,"w")
				self.writefile()
				return 1
		if valchk == 0:
			textmsg = ("Atencion, imposible guardar archivo, se encontraron %s errores" % self.totfall)
			msg = showerror(title = "Error", message = textmsg)
			msg2 = "Error, imposible guardar\n"
			self.txtout.insert("0.0", msg2)
			print msg2		
				
	def savedefault(self):
		"""Save default file, create a new file."""
		#check integrity--------------------------------------------
		valchk = self.checkseq()
		if valchk == 1:
			#escribe un fichero de salida
			self.cfgfile = "MDTA_default.cfg"		
			if self.cfgfile:
				self.txtfile =open(self.cfgfile,"w")
				self.writefile()	
				return 1
		if valchk == 0:
			textmsg = ("Atencion, imposible ejecutar aplicacion, se encontraron %s errores" % self.totfall)
			msg = showerror(title = "Error", message = textmsg)
			msg2 = "Error, imposible ejecutar\n"
			self.txtout.insert("0.0", msg2)
			print msg2			
	
	def writefile(self):
		"""Write the new created file with the GUI content."""
		#Corta las direcciones
		wkdir = os.getcwd()
		filename = self.dirin.get()
		finaldir = '/' + os.path.relpath(filename, wkdir)
		filename = self.dirout.get()
		finalout = '/' + os.path.relpath(filename, wkdir) + '/'
		#txtfile = open(self.cfgfile,"w")
		self.txtfile.writelines("VERSION %s\n" % self.currentv)
		self.txtfile.writelines("%s %s\n" % (self.txtout1, finaldir))
		self.txtfile.writelines("DIR_OUT %s\n" % finalout)
		#---------------------------------------------------------------
		self.txtfile.writelines("[%s]\n" % self.txtout2)
		if(self.newz.get() == 0):	
			self.txtfile.writelines("NEWZ 0\n")
			self.txtfile.writelines("%s 0\n" % self.txtout3)
			self.txtfile.writelines("%s 0\n" % self.txtout4)
			self.txtfile.writelines("%s /masknamefile.grd\n" % self.txtout5)
			self.txtfile.writelines("%s /xyzznamefile.txt\n" % self.txtout6)
		if(self.newz.get() == 1):	
			self.txtfile.writelines("NEWZ 1\n")	
			self.txtfile.writelines("%s %i\n" % (self.txtout3, self.nfase.get()))
			if self.nfase.get() == 1:
				
				filename = self.dirma.get()
				finalmask = '/' + os.path.relpath(filename, wkdir)
				
				self.txtfile.writelines("%s %i\n" % (self.txtout4, self.masktipe.get()))
				self.txtfile.writelines("%s %s\n" % (self.txtout5, finalmask))
				self.txtfile.writelines("%s /xyzznamefile.txt\n" % self.txtout6)
				
			if self.nfase.get() == 2:
				
				filename = self.nwz.get()
				finalxyz = '/' + os.path.relpath(filename, wkdir)
				
				self.txtfile.writelines("%s 0\n" % self.txtout4)
				self.txtfile.writelines("%s /masknamefile.grd\n" % self.txtout5)
				self.txtfile.writelines("%s %s\n" % (self.txtout6, finalxyz))
				
		#---------------------------------------------------------------		
		self.txtfile.writelines("[%s]\n" % self.txtout7)
		#self.txtfile.writelines("DEM_NAME %s\n" % self.demname.get())
		self.txtfile.writelines("%s %i\n" % (self.txtout8, self.demtipe.get()))
		self.txtfile.writelines("MAX_ZVAL %s\n" % self.maxval.get())
		self.txtfile.writelines("MIN_ZVAL %s\n" % self.minval.get())
		self.txtfile.writelines("NULL_VAL %i\n" % int(self.nullval.get()))
		if(self.recor.get() == 0):
			self.txtfile.writelines("%s 0\n" % self.txtout9)
			self.txtfile.writelines("X_MIN 0.0\n")
			self.txtfile.writelines("X_MAX 0.0\n")
			self.txtfile.writelines("Y_MIN 0.0\n")
			self.txtfile.writelines("Y_MAX 0.0\n")
			self.txtfile.writelines("RESX 0.0\n")
			self.txtfile.writelines("RESY 0.0\n")
		if(self.recor.get() == 1):
			self.txtfile.writelines("%s 1\n" % self.txtout9)	
			self.txtfile.writelines("X_MIN %s\n" % self.xmin.get())
			self.txtfile.writelines("X_MAX %s\n" % self.xmax.get())
			self.txtfile.writelines("Y_MIN %s\n" % self.ymin.get())
			self.txtfile.writelines("Y_MAX %s\n" % self.ymax.get())
			self.txtfile.writelines("RESX %s\n" % self.resx.get())  
			self.txtfile.writelines("RESY %s\n" % self.resy.get())
		#---------------------------------------------------------------
		self.txtfile.writelines("[%s]\n" % self.txtout10)
		valsik = int(self.sikvar.get())	
		self.txtfile.writelines("%s %i\n" % (self.txtout11, valsik))
		#aspcet ----------------------------
		valasp = int(self.aspvar.get())	
		self.txtfile.writelines("%s %i\n" % (self.txtout12, valasp))
		#slope	-----------------------------
		valslo = int(self.slopvar.get())	
		self.txtfile.writelines("%s %i\n" % (self.txtout13, valslo))		
		#---------------------------------------------------------------
		self.txtfile.writelines("[%s]\n" % self.txtout14)
		algo = int(self.tratip.get())
		if ( algo == 0):
			self.txtfile.writelines("%s 0\n" % self.txtout15)
			self.txtfile.writelines("DIST_MAX 0.0\n")
			self.txtfile.writelines("%s 0.0\n" % self.txtout16)
			self.txtfile.writelines("%s 0\n" % self.txtout17)
			self.txtfile.writelines("%s 0.0\n" % self.txtout18)
			self.txtfile.writelines("%s 0\n" % self.txtout19)
			self.txtfile.writelines("%s 0\n" % self.txtout22)       #<--------------
			self.txtfile.writelines("UTMZONE 0\n")
			self.txtfile.writelines("HEMIS 0\n")
			self.txtfile.writelines("[SEC_POINTS]\n")
			self.txtfile.writelines("N_CENTROS 0\n")
			self.txtfile.writelines("[X_Y_JERAR_RIO]\n")
		if ( algo > 0):
			self.txtfile.writelines("%s %s\n" % (self.txtout15, self.tratip.get()))
			self.txtfile.writelines("DIST_MAX %s\n" % self.dismax.get())
			self.txtfile.writelines("%s %s\n" % (self.txtout16, self.altmax.get()))
			if ( algo == 1 or algo == 2):
				self.txtfile.writelines("%s 0\n" % self.txtout17)
				self.txtfile.writelines("%s %s\n" % (self.txtout18, self.incre.get()))
				self.txtfile.writelines("%s %s\n" % (self.txtout23, self.force.get()))
				self.txtfile.writelines("%s 0\n" % self.txtout19)
			if ( algo == 3):		
				self.txtfile.writelines("%s 0\n" % self.txtout17)
				self.txtfile.writelines("%s %s\n" % (self.txtout18, self.incre.get()))
				self.txtfile.writelines("%s %s\n" % (self.txtout23, self.force.get()))
				self.txtfile.writelines("%s %s\n" % (self.txtout19, self.itera.get()))
			if ( algo == 4):
				self.txtfile.writelines("%s %s\n" % (self.txtout17, self.rad.get()))
				self.txtfile.writelines("%s 0.0\n" % self.txtout18)
				self.txtfile.writelines("%s 0\n" % self.txtout23)
				self.txtfile.writelines("%s 0\n" % self.txtout19)
			self.txtfile.writelines("WRITE_MOD %s\n" % self.mod.get())		
			self.txtfile.writelines("UTMZONE %s\n" % self.huso.get())
			self.txtfile.writelines("HEMIS %s\n" % self.hemis.get())
			self.txtfile.writelines("[SEC_POINTS]\n")
			self.txtfile.writelines("%s %s\n" % (self.txtout20, self.npt.get()))
			self.txtfile.writelines("[%s]\n" % self.txtout21)
			line = self.T.get("0.0", END)
			self.txtfile.writelines(line)		
		self.txtfile.close()
		msg2 = "End save file\n"		
		self.txtout.insert("0.0", msg2)
		print msg2
		return 1
		
	#LOAD FILE----------------------------------------------------------
	def loadfile(self):	
		#carga un fichero
		self.calload = 1
		"""Load an existed file in the GUI."""
		wkdir = os.getcwd()
		archivo = tkFileDialog.askopenfile(parent=self,mode='rb',title='Choose a cfg file')
		#archivo = tkFileDialog.askopenfilename(parent=self,title='Choose a cfg file')
		if archivo != None:
			self.txtout.delete("0.0",END)
			name = archivo.name
			cfgname = os.path.basename(name)
			textmsg = ("Leyendo archivo %s\n" % cfgname)
			self.txtout.insert("0.0", textmsg)
			self.cleangui()
			lines = archivo.readlines()
			lines = [line.rstrip('\n') for line in lines]
			#sec dir 0--------------------------------------------------
			val = lines[0].split(" ")
			vers = val[1]
			print self.currentv, vers
			if(self.currentv != vers):
				textmsg = ("Atencion, la version de archivo de configuracion no coincide")
				msg = showerror(title = "Error", message = textmsg)
				msg2 = ("Error en cfg\n")
				self.txtout.insert("0.0", msg2)
				print msg2
			else:	
				val = lines[1].split(" ")
				finaldir = wkdir + val[1]
				self.dirin.delete(0,END)
				self.dirin.insert(INSERT, finaldir)
				val = lines[2].split(" ")
				finalout = wkdir + val[1]
				self.dirout.delete(0,END)
				self.dirout.insert(INSERT, finalout)
				#sec modifica dem 3 ----------------------------------------
				val = lines[4].split(" ") #new dem
				dat = int(val[1])
				if(dat == 1):
					self.newz.set(1)
					self.update_chk()#----------CHECK	
					val = lines[5].split(" ") #fase
					dat2 = int(val[1])
					if(dat2 == 1):
						self.nfase.set(1)
					if(dat2 == 2):
						self.nfase.set(2)
				self.update_chk()#----------CHECK
				val = lines[6].split(" ")#tipo mascara
				dat = int(val[1])
				self.masktipe.set(dat)
				val = lines[7].split(" ")#mascara
				finalmask = wkdir + val[1]
				self.dirma.delete(0,END)
				self.dirma.insert(INSERT, finalmask)
				val = lines[8].split(" ")#xyzfile
				finalxyz = wkdir + val[1]
				self.nwz.delete(0,END)
				self.nwz.insert(INSERT, finalxyz)
				
				#val = lines[4].split(" ")
				#self.demname.delete(0,END)
				#self.demname.insert(INSERT, val[1])
				#secdem 9---------------------------------------------------
				val = lines[10].split(" ") #dem tipe
				dat = int(val[1])
				if(dat == 1):
					self.demtipe.set(1)
				if(dat == 2):
					self.demtipe.set(2)
				val = lines[11].split(" ") #max val
				self.maxval.delete(0,END)
				self.maxval.insert(INSERT, val[1])
				val = lines[12].split(" ") #min val
				self.minval.delete(0,END)
				self.minval.insert(INSERT, val[1])
				val = lines[13].split(" ")  #nul val
				self.nullval.delete(0,END)
				self.nullval.insert(INSERT, val[1])
				#recorte----------------------------
				val = lines[14].split(" ")  #recorte
				dat = int(val[1])
				if(dat == 1):	
					self.recor.set(1)
					self.update_chk()#----------CHECK
					val = lines[15].split(" ")  #xmin
					self.xmin.delete(0,END)
					self.xmin.insert(INSERT, val[1])
					val = lines[16].split(" ")  #xmax
					self.xmax.delete(0,END)
					self.xmax.insert(INSERT, val[1])
					val = lines[17].split(" ")  #ymin
					self.ymin.delete(0,END)
					self.ymin.insert(INSERT, val[1])
					val = lines[18].split(" ")  #ymax
					self.ymax.delete(0,END)
					self.ymax.insert(INSERT, val[1])
					val = lines[19].split(" ")  #resx
					self.resx.delete(0,END)
					self.resx.insert(INSERT, val[1])
					val = lines[20].split(" ")  #resy
					self.resy.delete(0,END)
					self.resy.insert(INSERT, val[1])
				#procesado 21-----------------------------------------------
				val = lines[22].split(" ")  #sink
				dat = int(val[1])
				if dat > 0 and dat < 3:
					self.sikvar.set(str(dat)) 
					self.update_chk()#----------CHECK
				if dat <= 0 or dat >= 3:	
					self.sikvar.set('0')
					self.update_chk()#----------CHECK	
						
				val = lines[23].split(" ") #asp
				dat = int(val[1])
				if dat > 0 and dat < 3:
					self.aspvar.set(str(dat)) 
					self.update_chk()#----------CHECK
				if dat <= 0 or dat >= 3:	
					self.aspvar.set('0')
					self.update_chk()#----------CHECK	
						
				val = lines[24].split(" ") #slope
				dat = int(val[1])
				if dat > 0 and dat < 10:
					self.slopvar.set(str(dat)) 
					self.update_chk()#----------CHECK
				if dat <= 0 or dat >= 10:	
					self.slopvar.set('0')
					self.update_chk()#----------CHECK

				#trayec 26--------------------------------------------------
				val = lines[26].split(" ")  #algoritmos
				self.tratip.set(val[1]) 
				alg = int(val[1])
				if(alg > 0):
					self.update_chk()#----------CHECK
					val = lines[27].split(" ")  #distmax
					self.dismax.delete(0,END)
					self.dismax.insert(INSERT, val[1])
					val = lines[28].split(" ")  #altmax
					self.altmax.delete(0,END)
					self.altmax.insert(INSERT, val[1])
					
					val = lines[33].split(" ")  #mod
					self.modval.delete(0,END)
					self.mod.set(val[1])
					
					val = lines[34].split(" ")  #huso
					self.huso.delete(0,END)
					self.huso.insert(INSERT, val[1])
					
					val = lines[35].split(" ")  #hemis
					self.hemisval.delete(0,END)
					self.hemis.set(val[1])                              #<--------------
					#self.hemisval.insert(INSERT, val[1])
					
					val = lines[37].split(" ")  #npt
					self.npt.delete(0,END)
					self.npt.insert(INSERT, val[1])
					npt = int(val[1])
					if(npt > 0):
						j=1;	
						for i in xrange(39, npt+39):  #puntos xy
							if(i == 39):
								self.T.delete("0.0",'end')
								self.T.insert("0.0", lines[i] + '\n')
							if(i > 39 and i < npt+39-1):	
								self.T.insert(str(j)+".0", lines[i] + '\n')
							if(i == npt+39-1):
								self.T.insert(str(j)+".0", lines[i])		
							j +=1;
					if(npt == 0):	
						self.npt.delete(0,END)
						self.npt.insert(INSERT, "0")		
					self.update_chk()#----------CHECK
					if(alg == 1 or alg == 2):	
						val = lines[30].split(" ")  #incremento
						self.incre.delete(0,END)
						self.incre.insert(INSERT, val[1])
						
						val = lines[31].split(" ")  #force
						self.forceval.delete(0,END)
						self.force.set(val[1])
						
						self.rad.delete(0,END)
						self.rad.insert(INSERT, "0.0")	
						self.itera.delete(0,END)
						self.itera.insert(INSERT, "0.0")	
					if(alg == 3):
						val = lines[30].split(" ")  #incremento
						self.incre.delete(0,END)
						self.incre.insert(INSERT, val[1])
							
						val = lines[31].split(" ")  #force
						self.forceval.delete(0,END)
						self.force.set(val[1])		
						
						val = lines[32].split(" ")  #itera
						self.itera.delete(0,END)
						self.itera.insert(INSERT, val[1])
						self.rad.delete(0,END)
						self.rad.insert(INSERT, "0.0")	
						#self.incre.delete(0,END)
						#self.incre.insert(INSERT, "0.0")
					if(alg == 4):
						val = lines[29].split(" ")  #radbus
						self.rad.delete(0,END)
						self.rad.insert(INSERT, val[1])	
						self.incre.delete(0,END)
						self.incre.insert(INSERT, "0.0")
						self.itera.delete(0,END)
						self.itera.insert(INSERT, "0.0")	
				#----------------------	
				archivo.close()	
				msg2 = "End load file\n"		
				self.txtout.insert("0.0", msg2)
				print msg2
				#check integrity--------------------------------------------
				msg2 = "Checking file\n"		
				self.txtout.insert("0.0", msg2)
				print msg2
				valchk = self.checkseq()
				if valchk == 0:
					textmsg = ("Atencion, se encontraron %s errores en la carga del archivo" % self.totfall)
					msg = showerror(title = "Error", message = textmsg)
					msg2 = ("%s Errores econtrados\n" % self.totfall)
					self.txtout.insert("0.0", msg2)
					print msg2
				if valchk == 1:	
					textmsg = "Archivo cargado correctamente\n"
					self.txtout.insert("0.0", textmsg)
					print msg2
	
	def loadnpt(self):	
		#carga un fichero
		"""Load an xyz file in the GUI."""
		wkdir = os.getcwd()
		archivo = tkFileDialog.askopenfile(parent=self,mode='rb',title='Choose a txt file')
		if archivo != None:
			self.txtout.delete("0.0",END)
			name = archivo.name
			cfgname = os.path.basename(name)
			textmsg = ("Leyendo archivo %s\n" % cfgname)
			self.txtout.insert("0.0", textmsg)
			#self.cleangui()
			self.T.delete("0.0",'end')
			lines = archivo.readlines()
			lines = [line.rstrip('\n') for line in lines]
			nlin = len(lines)
			#---
			val = lines[0].split()
			nval = len(val)
			if nval != 2:
				textmsg = "Atencion, el archivo debe tener 2 columnas - X Y - separadas por espacios"
				msg = showerror(title = "Error", message = textmsg)
				msg2 = "Error en estructura del archivo\n"
				self.txtout.insert("0.0", msg2)
				print msg2
				return 1
			if nval == 2:
				#check si hay cabecera
				self.valtx  = val[0]
				self.ndectype = 1
				resp1 = self.is_decimal()
				if(resp1 == 0):
					print "sin cabecera"
					ncab = 0
				if(resp1 == 1):	
					print "Con cabecera, saltamos linea"
					val = lines[1].split()
					ncab = 1

				self.valtx  = val[0]
				self.ndectype = 1
				resp1 = self.is_decimal()
				
				self.valtx  = val[1]
				self.ndectype = 1
				resp2 = self.is_decimal()
				
				#self.valtx  = val[2]
				#self.ndectype = 1
				#resp3 = self.is_entero()
				
				#self.valtx  = val[3]
				#self.ndectype = 1
				#resp4 = self.is_entero()
					
				#if resp1 == 0 and resp2 == 0 and resp3 == 0 and resp4 == 0:
				if resp1 == 0 and resp2 == 0:
					msg2 = "Estructura correcta, leyendo archivo\n"
					self.txtout.insert("0.0", msg2)
					print msg2
					j=0;
					for i in xrange(0, nlin):
						if(i == 0):
							if(ncab == 0):
								self.T.delete("0.0",'end')
								self.T.insert("1.0", lines[i] + '\n')
								j=2;
							if(ncab == 1):	
								j=1;		
						if(i > 0): 
							if(i < nlin-1):
								self.T.insert(str(j)+".0", lines[i] + '\n')	
							if(i == nlin-1):
								self.T.insert(str(j)+".0", lines[i])		
							j +=1;
					if ncab == 1:
						nlin = nlin - 1		
					self.npt.delete(0,END)
					self.npt.insert( INSERT, str(nlin) )	
					return 0
				else:
					textmsg = "Error en estructura, el archivo debe tener 2 columnas - X Y - separadas por espacios\n"
					msg = showerror(title = "Error", message = textmsg)
					self.txtout.insert("0.0", textmsg)
					print msg		
					return 1	
										
    #RUN MODEL---------------------------------------------------------- 
	def runmdt(self):
		savefile = self.savedefault()
		devuelve = 0
		if savefile != 1:	
			msg = "ATENCION, ocurrio un error, revise los datos e intentelo de nuevo\n"		
			self.txtout.insert("0.0", msg)
			print msg
		
		if savefile == 1:
			ejecutac = "./MDTanaliza MDTA_default.cfg"
			msg2 = "Inicio simulacion\n"		
			self.txtout.insert("0.0", msg2)
			print msg2
			os.system(ejecutac)
			devuelve = 1
			return 1
		if	devuelve == 1:
			msg2 = "Finalizacion simulacion\n"
			self.txtout.insert("0.0", msg2)
			print msg2
				
	#LANGUAJE GUI
	def lanGUI(self):
		#Language
		if (self.LAN == 1):
			self.txtbu0 = "Opciones"
			self.txtbu1 = "Cargar archivo cfg"
			self.txtbu2 = "Salvar archivo cfg"
			self.txtbu3 = "Ejecutar MDTanaliza"
			self.txtbu4 = "Ayuda"
			self.txtbu5 = "Salir"
			self.txtbu6 = "Change Lang Eng"
			
			self.txtgen1 = "DIRECTORIOS Y FICHEROS DE ENTRADA-SALIDA"
			self.txtgen2 = u"MODIFICACIÓN MDT"
			self.txtgen3 = u"PARÁMETROS MDT"
			self.txtgen4 = u"MORFOMETRÍA"
			self.txtgen5 = "TRAYECTORIA DE FLUJOS GRAVITACIONALES"
			self.txtgen6 = "PANEL INFORMATIVO"
			
			self.txtlab1 = "Formato MDT"
			self.txtlab2 = "Formato Mascara"
			self.txtlab3 = "Valor Z Max. del MDT (m)"
			self.txtlab4 = "Valor Z Min. del MDT (m)"
			self.txtlab5 = "Valor Nulo interno"
			self.txtlab6 = "X coordenada Min."
			self.txtlab7 = "X coordenada Max."
			self.txtlab8 = "Y coordenada Min."
			self.txtlab9 = "Y coordenada Max."
			self.txtlab10 = u"Resolución en X"
			self.txtlab11 = u"Resolución en Y"
			self.txtlab12 = u"Distancia máxima (m)"
			self.txtlab13 = u"Altura crítica (m)"
			self.txtlab14 = "Restric. Multiflow (%)"
			self.txtlab15 = "Incremento relleno (m)"
			self.txtlab16 = u"Iteraciones totales"
			self.txtlab17 = "Punto inicio totales"
			self.txtlab18 = "Tipo de algoritmo"
			self.txtlab18a = "Unidireccional Min.Top."
			self.txtlab18b = "Unidireccional Max.Pen."
			self.txtlab18c = "Montecarlo Random"
			self.txtlab18d = "Multidireccional Min.Top."
			self.txtlab19  = "Hemisferio"
			self.txtlab19a = "Hemisferio Norte"
			self.txtlab19b = "Hemisferio Sur"
			self.txtlab20  = "Modo Esc. Raster"
			self.txtlab21  = "Forzar Interacc."
			
			self.txtchk1 = "MDT modifica"
			self.txtchk2 = "Fase I"
			self.txtchk3 = "Fase II"
			self.txtchk4 = "Binario"
			self.txtchk5 = "ASCII  "
			self.txtchk6 = "Binario"
			self.txtchk7 = "ASCII  "
			self.txtchk8 = "Recorte - UTM (m)"
			self.txtchk9 = "Superficie sin salida"
			self.txtchk9a = u"S.S. Detección"
			self.txtchk9b = u"S.S. Modificación"
			self.txtchk10 = u"Pendiente-Orientación"
			self.txtchk10a = "P-O. Min.Top."
			self.txtchk10b = "P-O. Max.Pen."
			self.txtchk11 = "Pendiente-Gradiente"
			self.txtchk11a = "S-G Burrough and McDonell, 1998" #1
			self.txtchk11b = "S-G Fleming and Hoffer, 1979" #2
			self.txtchk11c = "S-G Unwin, 1981" #3
			self.txtchk11d = "S-G Sharpnack et al, 1969" #4
			self.txtchk11e = "S-G Horn, 1981"
			self.txtchk11f = "S-G Chu and Tsai, 1995"
			self.txtchk11g = "S-G Travis etal 1975, EPPL7, 1987"
			self.txtchk11h = "S-G Jones, 1998"
			self.txtchk11i = "S-G Wood, 1996"
			
			self.txtbut1 = "Archivo MDT"
			self.txtbut2 = "Directorio Salida"
			self.txtbut3 = u"Máscara"
			self.txtbut4 = "XYZ file"
			
			self.txtout1 = "MDT_IN"
			self.txtout2 = "SEC_MDT_MOD"
			self.txtout3 = "FASE"
			self.txtout4 = "FORMATO_MASCARA"
			self.txtout5 = "NOMBRE_MASCARA"
			self.txtout6 = "NOMBRE_XYZZ"
			self.txtout7 = "SEC_MDT"
			self.txtout8 = "FORMATO_MDT"
			self.txtout9 = "RECORTE"
			self.txtout10 = "MORFOMETRIA"
			self.txtout11 = "SIN_SALIDA"
			self.txtout12 = "PEND_ORIENT"
			self.txtout13 = "PEND_GRAD"
			self.txtout14 = "TRAY_FLUJOS"
			self.txtout15 = "ALGOR_TIPO"
			self.txtout16 = "CRIT_ALTURA"
			self.txtout17 = "REST_MULTIFLOW"
			self.txtout18 = "INCRE_RELLENO"
			self.txtout19 = "INTERACIONES"
			self.txtout20 = "PUNTOS_TOTALES"
			self.txtout21 = "PUNTOS_DATOS"
			self.txtout22 = "ESCRIB_MODO"
			self.txtout23 = "FORZAR_INTER"
			
			
		if (self.LAN == 2):
			self.txtbu0 = "Options"
			self.txtbu1 = "Load cfg file"
			self.txtbu2 = "Save cfg file"
			self.txtbu3 = "Run MDTanaliza"
			self.txtbu4 = "Help"
			self.txtbu5 = "Exit"
			self.txtbu6 = "Cambiar idioma Spa"
			
			self.txtgen1 = "INPUT-OUTPUT FILES AND DIRECTORIES"
			self.txtgen2 = "DEM MODIFICATION"
			self.txtgen3 = "DEM PARAMETERS"
			self.txtgen4 = "MORPHOMETRY"
			self.txtgen5 = "GRAVITY FLOW PAHTS"
			self.txtgen6 = "INFO PANEL"
			
			self.txtlab1 = "DEM Format"
			self.txtlab2 = "Mask Format"
			self.txtlab3 = "DEM Z Max. Value (m)"
			self.txtlab4 = "DEM Z Min. Value (m)"
			self.txtlab5 = "Internal Null value"
			self.txtlab6 = "X Min. coordenate"
			self.txtlab7 = "X Max. coordenate"
			self.txtlab8 = "Y Min. coordenate"
			self.txtlab9 = "Y Max. coordenate"
			self.txtlab10 = "X Resolution"
			self.txtlab11 = "Y Resolution"
			self.txtlab12 = "Maximum Distance (m)"
			self.txtlab13 = "Critical height (m)"
			self.txtlab14 = "Restric. Multiflow (%)"
			self.txtlab15 = "Fill increase (m)"
			self.txtlab16 = "Total Iterations"
			self.txtlab17 = "Total init. points"
			self.txtlab18 = "Algorithm type"
			self.txtlab18a = "Single Path LHM"
			self.txtlab18b = "Single Path SSM"
			self.txtlab18c = "Montecarlo Random"
			self.txtlab18d = "Multiple Path. LHM"
			self.txtlab19  = "Hemisphere"
			self.txtlab19a = "Northern Hemisphere"
			self.txtlab19b = "Southern Hemisphere"
			self.txtlab20  = "W. Raster Mode"
			self.txtlab21  = "Force Interacc."
			
			self.txtchk1 = "Modified DEM"
			self.txtchk2 = "Step I"
			self.txtchk3 = "Step II"
			self.txtchk4 = "Binary"
			self.txtchk5 = "ASCII  "
			self.txtchk6 = "Binary"
			self.txtchk7 = "ASCII  "
			self.txtchk8 = "Clip - UTM (m)"
			self.txtchk9 = "Surface Depression"
			self.txtchk9a = "S.D. Detection"
			self.txtchk9b = "S.D. Modification"
			
			self.txtchk10 = "Slope-Aspect"
			self.txtchk10a = "S-A. LHM"
			self.txtchk10b = "S-A. SSM"
			
			self.txtchk11 = "Slope-Gradient"
			self.txtchk11a = "S-G Burrough and McDonell, 1998" #1
			self.txtchk11b = "S-G Fleming and Hoffer, 1979" #2
			self.txtchk11c = "S-G Unwin, 1981" #3
			self.txtchk11d = "S-G Sharpnack et al, 1969" #4
			self.txtchk11e = "S-G Horn, 1981"
			self.txtchk11f = "S-G Chu and Tsai, 1995"
			self.txtchk11g = "S-G Travis etal 1975, EPPL7, 1987"
			self.txtchk11h = "S-G Jones, 1998"
			self.txtchk11i = "S-G Wood, 1996"
			
			self.txtbut1 = "DEM File"
			self.txtbut2 = "Output directory"
			self.txtbut3 = "Raster Mask"
			self.txtbut4 = "XYZ file"
			
			self.txtout1 = "DEM_IN"
			self.txtout2 = "SEC_MOD_DEM"
			self.txtout3 = "PHASE"
			self.txtout4 = "MASK_FORMAT"
			self.txtout5 = "MASK_FILENAME"
			self.txtout6 = "XYZZ_FILENAME"
			self.txtout7 = "SEC_DEM"
			self.txtout8 = "DEM_FORMAT"
			self.txtout9 = "CLIP"
			self.txtout10 = "MORPHOMETRY"
			self.txtout11 = "SURF_DEPRESSION"
			self.txtout12 = "SLOPE_ASPECT"
			self.txtout13 = "SLOPE_GRAD"
			self.txtout14 = "FLOW_PATH"
			self.txtout15 = "ALGOR_TYPE"
			self.txtout16 = "CRIT_HEIGHT"
			self.txtout17 = "REST_MULTIFLOW"
			self.txtout18 = "FILL_INCRE"
			self.txtout19 = "ITERATIONS"
			self.txtout20 = "TOTAL_POINTS"
			self.txtout21 = "POINT_DATA"
			self.txtout22 = "WRITE_MOD"
			self.txtout23 = "FORCE_INTER"
	
	#GRAFIC GUI---------------------------------------------------------
	def initUI(self):
		self.ndectype = 0
		self.ntipe = 0
		self.parent.title("MDTanaliza v. 2.1-2018-04-06")
		self.currentv = '1.2'
		self.pack(fill=BOTH, expand=1)
		self.lanGUI()
		
		#***************************************************************
		#MENU BAR-------------------------------------------------------
		
		self.menubar = Menu(self)
		menu = Menu(self.menubar, tearoff=0)
		self.menubar.add_cascade(label=self.txtbu0, menu=menu)
		if (self.LAN == 2):menu.add_command(label=self.txtbu6, comman=self.chalangen2sp)
		if (self.LAN == 1):menu.add_command(label=self.txtbu6, comman=self.chalangsp2en)
		menu.add_command(label=self.txtbu1, comman=self.loadfile)
		menu.add_command(label=self.txtbu2, comman=self.savefile)
		menu.add_command(label=self.txtbu3, comman=self.runmdt)
		menu.add_command(label=self.txtbu4, comman=self.openhelp)
		menu.add_command(label=self.txtbu5, command=self.quit)

		try:
			self.master.config(menu=self.menubar)
		except AttributeError:
			# master is a toplevel window (Python 1.4/Tkinter 1.63)
			self.master.tk.call(master, "config", "-menu", self.menubar)			
		
        #***************************************************************
        #SECCION DIRECTORIOS directorios--------------------------------
        #labels
		labdirgen = Label(self, text=self.txtgen1, fg = "red")
		labdirgen.grid(row=0, columnspan=4, sticky=W+E)
		labnwz = Label(self, text=self.txtgen2, fg = "red")
		labnwz.grid(row=3, columnspan=4, sticky=W+E)
		#campos
		self.dirin   = Entry(self)
		self.dirin.grid(row=1, column=1, columnspan=2, sticky=W+E)
		self.dirout  = Entry(self)
		self.dirout.grid(row=2, column=1, columnspan=2, sticky=W+E)
		#---
		self.newz = BooleanVar()
		self.newbut = Checkbutton(self, text=self.txtchk1, variable=self.newz, command = self.update_chk)
		self.newbut.grid(row=3, column=3)
		self.dirma   = Entry(self, state=DISABLED)
		self.dirma.grid(row=4, column=1, columnspan=2, sticky=W+E)
		#botones
		self.demButton = Button(self, text=self.txtbut1, command=self.dtmfile, width=20)
		self.demButton.grid(row=1)
		self.outButton = Button(self, text=self.txtbut2, command=self.direout, width=20)
		self.outButton.grid(row=2)
		#---
		self.maskButton = Button(self, text=self.txtbut3, width=20, state=DISABLED, command=self.maskfile)
		self.maskButton.grid(row=4)
		self.newButton = Button(self, text=self.txtbut4, command=self.newopenfile, state=DISABLED, width=20)
		self.newButton.grid(row=5)
		self.nwz     = Entry(self, justify=LEFT, width=25, state=DISABLED)
		self.nwz.grid(row=5, column=1, columnspan=2, sticky=W+E)
		self.nfase = IntVar()
		self.fase1 = Radiobutton(self, 
			text=self.txtchk2,
			#padx = 20, 
			justify=LEFT,
			variable=self.nfase, 
			width=25,
			value=1,
			state=DISABLED, command = self.update_chk)
		self.fase1.grid(row=4, column=3)
		self.fase2 = Radiobutton(self, 
			text=self.txtchk3,
			justify=LEFT,
			#padx = 20, 
			variable=self.nfase,
			width=25, 
			value=2,
			state=DISABLED, command = self.update_chk)
		self.fase2.grid(row=5, column=3)
		
		logo = PhotoImage(file="MDTanaliza_ico.gif",)
		labimg = Label(self, image=logo, relief=RIDGE, height=65, width=65)
		labimg.logo = logo
		labimg.grid(row=0, column=3, rowspan=3)
		
		#***************************************************************
		#SECCION DATOS DTM----------------------------------------------
		#labels EMPIEZA EN 6 - 24
		labdemgen = Label(self, text=self.txtgen3, fg = "red", width=25).grid(row=6, columnspan=2, sticky=W+E)
		self.labdtipe = Label(self, justify=LEFT, text=self.txtlab1)
		self.labdtipe.grid(row=7)
		self.labmasktipe = Label(self, justify=LEFT, text=self.txtlab2)
		self.labmasktipe.grid(row=7, column=1)
		self.labdmax  = Label(self, justify=LEFT, text=self.txtlab3, relief=RIDGE, width=25)
		self.labdmax.grid(row=10)
		self.labdmin  = Label(self, justify=LEFT, text=self.txtlab4, relief=RIDGE, width=25)
		self.labdmin.grid(row=11)
		self.labdnull = Label(self, justify=LEFT, text=self.txtlab5, relief=RIDGE, width=25)
		self.labdnull.grid(row=12)
		self.labdxmin = Label(self, justify=LEFT, text=self.txtlab6, relief=RIDGE, width=25)
		self.labdxmin.grid(row=14)
		self.labdxmax = Label(self, justify=LEFT, text=self.txtlab7, relief=RIDGE, width=25)
		self.labdxmax.grid(row=14, column=1)
		self.labdymin = Label(self, justify=LEFT, text=self.txtlab8, relief=RIDGE, width=25)
		self.labdymin.grid(row=16)
		self.labdymax = Label(self, justify=LEFT, text=self.txtlab9, relief=RIDGE, width=25)
		self.labdymax.grid(row=16, column=1)
		self.labdrex = Label(self, justify=LEFT, text=self.txtlab10, relief=RIDGE, width=25)
		self.labdrex.grid(row=18)
		self.labdrey = Label(self, justify=LEFT, text=self.txtlab11, relief=RIDGE, width=25)
		self.labdrey.grid(row=18, column=1)
		#campos
		self.demtipe = IntVar()
		self.demradb = Radiobutton(self, text=self.txtchk4, justify=LEFT, variable=self.demtipe, command = self.checkdirectories, width=25, value=1)
		self.demradb.grid(row=8)
		self.demrada = Radiobutton(self, text=self.txtchk5, justify=LEFT, variable=self.demtipe, command = self.checkdirectories, width=25, value=2)
		self.demrada.grid(row=9)
		
		self.masktipe = IntVar()
		self.maskradb = Radiobutton(self, text=self.txtchk6, justify=LEFT, variable=self.masktipe, command = self.checkdirectories, state=DISABLED, width=25, value=1)
		self.maskradb.grid(row=8, column=1)
		self.maskrada = Radiobutton(self, text=self.txtchk7, justify=LEFT, variable=self.masktipe, command = self.checkdirectories, state=DISABLED, width=25, value=2)
		self.maskrada.grid(row=9, column=1)
		
		self.maxval = Entry(self, justify=RIGHT, width=25)
		self.maxval.grid(row=10, column=1)
		self.minval = Entry(self, justify=RIGHT, width=25)
		self.minval.grid(row=11, column=1)
		self.nullval = Entry(self, justify=RIGHT, width=25)
		self.nullval.grid(row=12, column=1)
		self.nullval.insert(INSERT, "-9999")
		self.recor = BooleanVar()
		self.recorbut = Checkbutton(self, text=self.txtchk8, variable=self.recor,  command = self.update_chk)
		self.recorbut.grid(row=13, sticky=W)
		self.xmin = Entry(self, justify=RIGHT, width=25, state=DISABLED)
		self.xmax = Entry(self, justify=RIGHT, width=25, state=DISABLED)
		self.ymin = Entry(self, justify=RIGHT, width=25, state=DISABLED)
		self.ymax = Entry(self, justify=RIGHT, width=25, state=DISABLED)
		self.resx = Entry(self, justify=RIGHT, width=25, state=DISABLED)
		self.resy = Entry(self, justify=RIGHT, width=25, state=DISABLED)
		self.xmin.grid(row=15)
		self.xmax.grid(row=15, column=1)
		self.ymin.grid(row=17)
		self.ymax.grid(row=17, column=1)
		self.resx.grid(row=19)
		self.resy.grid(row=19, column=1)		
		
		#***************************************************************
		#SECCION OPCIONES PROCESADO DTM---------------------------------
		#labels EMPIEZA EN 20
		labprogen = Label(self, text=self.txtgen4, fg = "red", width=25).grid(row=20, columnspan=2, sticky=W+E)
		#campos
		#SINK
		self.siklab  = Label(self, justify=LEFT, text=self.txtchk9, relief=RIDGE, width=27)
		self.siklab.grid(row=21)
		self.sikvar = StringVar()
		self.sikval = Spinbox(self, from_=0, to=2, textvariable=self.sikvar, command = self.update_chk, state="readonly")
		self.sikval.grid(row=21, column=1, sticky=W+E)
		#ASPECT
		self.asplab  = Label(self, justify=LEFT, text=self.txtchk10, relief=RIDGE, width=27)
		self.asplab.grid(row=22)
		self.aspvar = StringVar()
		self.aspval = Spinbox(self, from_=0, to=2, textvariable=self.aspvar, command = self.update_chk, state="readonly")
		self.aspval.grid(row=22, column=1, sticky=W+E)
		#SLOPE
		self.sloplab  = Label(self, justify=LEFT, text=self.txtchk11, relief=RIDGE, width=27)
		self.sloplab.grid(row=23) 
		self.slopvar = StringVar()
		self.slopval = Spinbox(self, from_=0, to=9, textvariable=self.slopvar, command = self.update_chk, state="readonly")
		self.slopval.grid(row=23, column=1, sticky=W+E) 
    
		#***************************************************************
		#SECCION TRAYECTORIAS-------------------------------------------
		#labels EMPIEZA EN 6
		self.labdminlabtragen = Label(self, text=self.txtgen5, fg = "red", width=25)
		self.labdminlabtragen.grid(row=6, column=2, columnspan=2, sticky=W+E)
		self.labalgtyp = Label(self,text=self.txtlab18, background='yellow', foreground="blue", relief=RIDGE, width=25)
		self.labalgtyp.grid(row=7, column=2)
		self.labdminlabdmax = Label(self, justify=LEFT, text=self.txtlab12, relief=RIDGE, width=25)
		self.labdminlabdmax.grid(row=8,column=2)
		self.labdminlabalt = Label(self, justify=LEFT, text=self.txtlab13, relief=RIDGE, width=25)
		self.labdminlabalt.grid(row=9,column=2)
		self.labdminlabrad = Label(self, justify=LEFT, text=self.txtlab14, relief=RIDGE, width=25)
		self.labdminlabrad.grid(row=10,column=2)
		self.labdminlabinc = Label(self, justify=LEFT, text=self.txtlab15, relief=RIDGE, width=25)
		self.labdminlabinc.grid(row=11,column=2)
		self.labdminlabfor = Label(self, justify=LEFT, text=self.txtlab21, relief=RIDGE, width=25) #<---------------
		self.labdminlabfor.grid(row=12,column=2)
		self.labdminlabite = Label(self, justify=LEFT, text=self.txtlab16, relief=RIDGE, width=25)
		self.labdminlabite.grid(row=13,column=2)
		self.labdminlabmod = Label(self, justify=LEFT, text=self.txtlab20, relief=RIDGE, width=25)   #<---------------
		self.labdminlabmod.grid(row=14,column=2)
		
		self.labdminlabnpt = Label(self, justify=LEFT, text=self.txtlab17, relief=RIDGE, width=25)
		self.labdminlabnpt.grid(row=15,column=2)
		#campos	
		self.tratip  = StringVar()
		self.trayval = Spinbox(self, from_=0, to=4, textvariable=self.tratip, width=20, command = self.update_chk, state="readonly")
		self.dismax  = Entry(self, justify=RIGHT, width=20, state=DISABLED)
		self.altmax  = Entry(self, justify=RIGHT, width=20, state=DISABLED)
		self.rad     = Entry(self, justify=RIGHT, width=20, state=DISABLED)
		self.incre   = Entry(self, justify=RIGHT, width=20, state=DISABLED)
		
		#self.force   = Entry(self, justify=RIGHT, width=20, state=DISABLED)           #<---------------
		self.force  = StringVar()
		self.forceval = Spinbox(self, from_=0, to=1, textvariable=self.force, width=20, command = self.update_chk, state=DISABLED)
		
		self.itera   = Entry(self, justify=RIGHT, width=20, state=DISABLED)
		
		#self.mod     = Entry(self, justify=RIGHT, width=20, state=DISABLED)            #<---------------
		self.mod  = StringVar()
		self.modval = Spinbox(self, from_=0, to=1, textvariable=self.mod, width=20, command = self.update_chk, state=DISABLED)
		
		self.npt = Entry(self, justify=RIGHT, width=20)
		self.npt.delete(0,END)
		self.npt.insert( INSERT, "0" )
		self.trayval.grid(row=7,column=3)
		self.dismax.grid(row=8, column=3)
		self.altmax.grid(row=9, column=3)
		self.rad.grid(row=10, column=3)
		self.incre.grid(row=11, column=3)
		self.forceval.grid(row=12, column=3)                                            #<---------------
		self.itera.grid(row=13, column=3)
		self.modval.grid(row=14, column=3)                                                #<---------------
		self.npt.grid(row=15, column=3)  #de 13 a 14
		
		
		#LOAD XY VENTS
		self.S = Scrollbar(self)
		self.T = Text(self, height=4, width=27, state=DISABLED)
		self.S.grid(row=15, column=3, rowspan=5, sticky=N+S+E)
		#self.T.grid(row=14, column=3, rowspan=5, columnspan=1, pady=3, sticky=W)
		self.T.grid(row=15, column=3, rowspan=5)
		self.S.config(command=self.T.yview)
		self.T.config(yscrollcommand=self.S.set)
		self.loadnpt = Button(self, text="Load xy file", width=20, comman=self.loadnpt, state=DISABLED)
		self.loadnpt.grid(row=16, column=2,)
		#--
		self.labhuso = Label(self, justify=LEFT, text="UTM ZONE", relief=RIDGE, width=20)
		self.labhuso.grid(row=20,column=3)
		self.huso = Entry(self, justify=RIGHT, width=20, state=DISABLED)
		self.huso.grid(row=21, column=3,)
		self.labhemis = Label(self, justify=LEFT, text=self.txtlab19, relief=RIDGE, width=20)
		self.labhemis.grid(row=22,column=3)
		#self.hemis = Entry(self, justify=RIGHT, width=20, state=DISABLED)
		#self.hemis.grid(row=23, column=3,)
		
		self.hemis  = StringVar()
		self.hemisval = Spinbox(self, from_=0, to=1, textvariable=self.hemis, width=20, command = self.update_chk, state=DISABLED)
		self.hemisval.grid(row=23, column=3,)
		
		#Salida de mensajes---------------------------------------------
		self.labpanel = Label(self, text=self.txtgen6, fg = "red", width=25)
		self.labpanel.grid(row=17,column=2)
		
		self.scr = Scrollbar(self)
		self.txtout = Text(self, height=8, width=26)
		self.scr.grid(row=18, column=2, rowspan=6, sticky=N+S+E)
		self.txtout.grid(row=18, column=2, rowspan=6, sticky=W)
		self.scr.config(command=self.txtout.yview)
		self.txtout.config(yscrollcommand=self.scr.set)
		


def main():
  
    root = Tk()
    root.geometry("920x550+300+300")
    #root.state('zoomed')
    app = Example(root)
    root.mainloop()  


if __name__ == '__main__':
    main()
