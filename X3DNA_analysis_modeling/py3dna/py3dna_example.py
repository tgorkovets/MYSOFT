#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Простой пример использования py3dna
'''
from py3dna import py3dna
# Класс написан так, что может работать из VMD
# и всемсто PDB принимать atomsel (VMD_atomsel=atomsel), может быть удобно
# для расчета энергии деформации ДНК по трактории
# Инициализируем  объект из PDB
dna=py3dna('1KX5.pdb',tempdir='temp/'
    ,path='/home/armeev/Software/Source/x3dna-v2.1')
    
# или в VMD
#dnaSel=atomsel('nucleic')
#dna=py3dna(VMD_atomsel=dnaSel,tempdir='temp/'
#    ,path='/home/armeev/Software/Source/x3dna-v2.1')

# Включаем участки ДНК в процедуру минимизации
# Здесь выбираем участки с 1 по 40 и 107 по 147 нп включительно
# шаги между ними будут изменяться при минимизации
dna.set_movable_bp([[1,40],[107,147]])

#Выбираем расстояния, которые будут ограничиваться при минимизации
# Здесь мы ограничиваем расстояние между 1 и 147 нп 80 А
# и расстояние между 10 и 137 нп 75 А
dna.set_pairs_list_and_dist([[1,147],[10,437]],[80,75])

#Минимизируем конформацию без учета расстояний
frame,result=dna.run_minimize(usepairs=False,verbose=True)
# пишем конформацию в PDB
# путь должен быть абсолютным, или файл будет записан в TEMP
dna.frame_to_pdb(frame,'/home/armeev/result.pdb')

#Продолжаем минимизацию, включив ограничнеия по парам
#Ограничив количество итераций алгоритма минимизации
frame,result=dna.run_minimize(frame=frame,usepairs=True,verbose=True
    ,maxiter=20,maxfev=20)
dna.frame_to_pdb(frame,'/home/armeev/result_dist.pdb')
