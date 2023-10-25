import sys
#Usage: python3 path/to/csv [-s save.csv]


args = sys.argv[1:]
if '-s' in args:
    saveFn = args.pop(args.index('-s')+1)
else:
    saveFn = None
fns = [fn for fn in args if not fn.startswith('-')]


standardAmino = {'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR'}
standardDNA = {'DA', 'DC', 'DG', 'DT', 'DI'}
standardRNA = {'A', 'C', 'G','U', 'I'}
standardHOH = {'HOH'}

fnallRes = set()
fncountDict = dict() #dict of dicts {fn: resDict}

for fn in fns:
    resDict = dict() #dict {res:set of ids}
    pdbFile = open(fn, 'r')
    lines = [line for line in pdbFile.readlines() if line.startswith(('ATOM', 'HETATM'))]
    for line in lines:
        res=line[17:20].strip()
        resid=line[23:26].strip()
        if res in resDict.keys():
            resDict[res].add(resid)
        else:
            resDict[res] = {resid}
    allRes = set(resDict.keys())
    fnallRes.update(allRes)
    fncountDict[fn] = resDict
    nonstandardRes = allRes.difference(standardAmino, standardDNA, standardRNA, standardHOH)
    aminoRes = allRes.intersection(standardAmino)
    DNARes = allRes.intersection(standardDNA)
    RNARes = allRes.intersection(standardRNA)
    HOHRes = allRes.intersection(standardHOH)
    print(fn)
    if len(aminoRes) > 0: 
        print('--Amino Acids--')
        for res in sorted(list(aminoRes)):
            print('{}:{}'.format(res, len(resDict[res])))
    if len(DNARes) > 0: 
        print('--DNA--')
        for res in sorted(list(DNARes)):
            print('{}:{}'.format(res, len(resDict[res])))
    if len(RNARes) > 0: 
        print('--RNA--')
        for res in sorted(list(RNARes)):
            print('{}:{}'.format(res, len(resDict[res])))
    if len(nonstandardRes) > 0: 
        print('--Nonstandard--')
        for res in sorted(list(nonstandardRes)):
            print('{}:{}'.format(res, len(resDict[res])))
    if len(HOHRes) > 0: 
        print('--Water--')
        for res in HOHRes:
            print('{}:{}'.format(res, len(resDict[res])))
    
if not isinstance(saveFn, type(None)):
    saveFile = open(saveFn, 'w')
    fnList = sorted(list(fncountDict.keys()))
    ln = ''
    for fn in fnList: ln += ',{}'.format(fn)
    ln+='\n'
    saveFile.write(ln)

    nonstandardRes = sorted(list(fnallRes.difference(standardAmino, standardDNA, standardRNA, standardHOH)))
    aminoRes =       sorted(list(fnallRes.intersection(standardAmino)))
    DNARes =         sorted(list(fnallRes.intersection(standardDNA)))
    RNARes =         sorted(list(fnallRes.intersection(standardRNA)))
    HOHRes =         sorted(list(fnallRes.intersection(standardHOH)))
    resList = aminoRes + DNARes + RNARes + nonstandardRes + HOHRes
    for res in resList:
        ln = '{}'.format(res)
        for fn in fnList:
            resDict = fncountDict[fn]
            if res in resDict:
                ln += ',{}'.format(len(resDict[res]))
            else:
                ln +=',0'
        ln += '\n'
        saveFile.write(ln)
    saveFile.close()

