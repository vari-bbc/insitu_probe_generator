def df_select (name,fullseq,amplifier,pause,choose,polyAT,polyCG,BlastProbes,db,dropout,show,report,maxprobe,numbr):
    from Bio.Seq import Seq
    from Bio.Blast.Applications import NcbiblastnCommandline as bn
    import io
    import numpy as np
    import pandas as pd
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.max_rows',5000)
    pd.set_option('display.width', 80)
    from datetime import date
    from IPython.display import display, Markdown
    from ipydatagrid import DataGrid, TextRenderer, Expr

    __author__ = "Ryan W Null - ORCID_0000-0002-3830-4152"
    __copyright__ = "Copyright 2019-2021  The Ozpolat Lab,  https://bduyguozpolat.org/ "
    __credits__ = ['M. Desmond Ramirez - ORCID_0000-0003-3873-6999','Dennis Sun - ORCID_0000-0003-1000-7276','B. Duygu Ozpolat - ORCID_0000-0002-1900-965X']
    __source__ = "https://github.com/rwnull/insitu_probe_generator"
    __license__ = "GPL 3.0"
    __version__ = "2021_0.3.2"
    __DOI__ = "https://doi.org/10.5281/zenodo.3871970"

    def amp(ampl): 
        if ampl == "B1":
            upspc= "aa"
            dnspc= "ta"
            up = "GAGGAGGGCAGCAAACGG"
            dn = "GAAGAGTCTTCCTTTACG"
        elif ampl == "B2":
            upspc= "aa"
            dnspc= "aa"
            up = "CCTCGTAAATCCTCATCA"
            dn = "ATCATCCAGTAAACCGCC"
        elif ampl == "B3":
            upspc= "tt"
            dnspc= "tt"
            up = "GTCCCTGCCTCTATATCT"
            dn = "CCACTCAACTTTAACCCG"
        elif ampl == "B4":
            upspc= "aa"
            dnspc= "at"
            up = "CCTCAACCTACCTCCAAC"
            dn = "TCTCACCATATTCGCTTC"
        elif ampl == "B5":
            upspc= "aa"
            dnspc= "aa"
            up = "CTCACTCCCAATCTCTAT"
            dn = "CTACCCTACAAATCCAAT"
        elif ampl == "B7":
            upspc= "ww"
            dnspc= "ww"
            up = "CTTCAACCTCCACCTACC"
            dn = "TCCAATCCCTACCCTCAC"
        elif ampl == "B9":
            upspc= "ww"
            dnspc= "ww"
            up = "CACGTATCTACTCCACTC"
            dn = "TCAGCACACTCCCAACCC"
        elif ampl == "B10":
            upspc= "ww"
            dnspc= "ww"
            up = "CCTCAAGATACTCCTCTA"
            dn = "CCTACTCGACTACCCTAG"
        elif ampl == "B11":
            upspc= "ww"
            dnspc= "ww"
            up = "CGCTTAGATATCACTCCT"
            dn = "ACGTCGACCACACTCATC"
        elif ampl == "B13":
            upspc= "ww"
            dnspc= "ww"
            up = "AGGTAACGCCTTCCTGCT"
            dn = "TTATGCTCAACATACAAC"
        elif ampl == "B14":
            upspc= "ww"
            dnspc= "ww"
            up = "AATGTCAATAGCGAGCGA"
            dn = "CCCTATATTTCTGCACAG"
        elif ampl == "B15":
            upspc= "ww"
            dnspc= "ww"
            up = "CAGATTAACACACCACAA"
            dn = "GGTATCTCGAACACTCTC"
        elif ampl == "B17":
            upspc= "ww"
            dnspc= "ww"
            up = "CGATTGTTTGTTGTGGAC"
            dn = "GCATGCTAATCGGATGAG"
        else:
            print ("Please try again")
        return([upspc,dnspc,up,dn])

    def max33(maxprobe, seqs, numbr):
        """
        Filter `seqs` down to at most `numbr` items using the same spacing logic
        as your original function. Returns a dictionary with:
        - reduced:   the filtered list of seqs
        - keep_idx:  ndarray of original positions kept (use for df.iloc)
        - removed_idx: ndarray of original positions removed
        - mask:      boolean mask of length len(seqs) (True = kept)
        - message:   optional info string (e.g., when no action taken)
        """
        seqs = list(seqs)
        n = len(seqs)
        num = int(numbr)

        # Default: no change
        result = {
            "reduced": list(seqs),
            "keep_idx": np.arange(n, dtype=int),
            "removed_idx": np.array([], dtype=int),
            "mask": np.ones(n, dtype=bool),
        }

        # Only reduce when requested
        if maxprobe == "Yes" and num < n:
            keep = 33 if num == 0 else num
            keep = max(0, min(keep, n))

            # How many items will be skipped (removed)
            skip = n - keep
            if keep == 0:
                # Keep nothing
                mask = np.zeros(n, dtype=bool)
                result.update({
                    "reduced": [],
                    "keep_idx": np.array([], dtype=int),
                    "removed_idx": np.arange(n, dtype=int),
                    "mask": mask,
                })
                return result

            # Distribute zeros (skips) between ones (kept) exactly like your code
            zeros_per_one = skip // keep
            extra_zeros = skip - keep * zeros_per_one  # first `extra_zeros` gaps get one extra zero

            kept_positions = []
            pos = 0
            used_extra = 0
            for _ in range(keep):
                # place a "1" (keep) at current position
                kept_positions.append(pos)
                pos += 1
                # optional extra zero
                if used_extra < extra_zeros:
                    pos += 1
                    used_extra += 1
                # the regular zeros
                pos += zeros_per_one

            kept_positions = np.asarray(kept_positions, dtype=int)

            mask = np.zeros(n, dtype=bool)
            mask[kept_positions] = True
            removed_idx = np.where(~mask)[0]

            reduced = [seqs[i] for i in kept_positions]

            result.update({
                "reduced": reduced,
                "keep_idx": kept_positions,
                "removed_idx": removed_idx,
                "mask": mask,
            })
            return result
        return result



    # def output(cdna,g,fullseq,count,amplifier,name,pause,seqs,df_grid):

    #     if int(count) > 0:
    #         print()
    #         print()
    #         print("Figure Layout of Probe Sequences:")
    #         print("")
    #         print(str(amplifier+"_"+str(name)+"_PP"+str(count)+"_Dla"+str(pause)))
    #         grid.auto_fit_columns = True
    #         display(grid)
    #         print()
    #         print("Please select which probes you would like to keep")
    #         print()
    #         print()
    #         print()
    #         print()
    #         print("This is the in-place localization of the probe pairs along the full-length sense cDNA.")
    #         print()
    #         print(">"+name+" Sense Strand")
    #         print(g)
    #         print()
    #         print()
    #         print()
    #         print()
    #         print("Anti-sense sequence used to create probes:")
    #         print()
    #         print(">"+name+" Anti-Sense Strand")
    #         print(fullseq)    
    #         return()




    ### Printing out header
    
    display(Markdown("## HCR3.0 Probe Maker Output"))
    print("version" + __version__)
    print(__DOI__)
    print()
    print("Written by "+__author__) 
    print(" with "+__credits__[0])
    print(__credits__[1])
    print(' and '+__credits__[2])
    print()
    print(__copyright__)
    print(" with licensing provided under "+__license__)
    print()
    print("For more information visit: ")
    print(" "+__source__)
    print()
    print(date.today())
    print()

    display(Markdown("---"))    

    ### Init vars

    name=str(name)
    

    fullseq = Seq(fullseq)
    fullseq = fullseq.reverse_complement()
    fullseq = str(fullseq)
    cdna = len(fullseq)
    pause = int(pause)    
    
    
    amplifier=str((amplifier).upper())
    test=amp(amplifier)
    uspc=test[0]
    dspc=test[1]
    upinit=test[2]
    dninit=test[3]

    

    hpA = "A"*(polyAT+1)
    hpT = "T"*(polyAT+1)
    hpC = "C"*(polyCG+1)
    hpG = "G"*(polyCG+1)

    position = cdna-pause
    start = np.arange(0,cdna-52,1)
    end = np.arange(52,cdna,1)
    table = np.vstack([start,end])

    ### what are these?
    seqs={}
    pos=[]

    ### walk full sequence by a 52 bp window, ignore pos with homopolymers
    a=0

    while a < (position-52):
        if ((str(fullseq[table[0][a]:table[1][a]])).find(hpA) + (str(fullseq[table[0][a]:table[1][a]])).find(hpT) + (str(fullseq[table[0][a]:table[1][a]])).find(hpC) + (str(fullseq[table[0][a]:table[1][a]])).find(hpG) > -4):
            a += 1
        else:
            pos.append([table[0][a],table[1][a]])
            a += 1

    ### Creating the first trace through the sequence looking for max number of probe sequences 
    
    a = 0
    newlist = []
    newlista = []
    newlista2 = []
    newlistb = []
    newlistb2 = []
    strt=pos[0][0]
    stp=pos[0][1]
    newlista2.append([cdna-strt,cdna-stp])
    newlista.append([strt,stp])
    while a < len(pos):    
        if pos[a][0] > (stp + 2):
            strt = pos[a][0]
            stp  = pos[a][1]
            newlista2.append([cdna-strt,cdna-stp])
            newlista.append([strt,stp])
            a+=1
        else :
            a+=1
    lists = {}
    listz = {}
    listz[0] = newlista2
    lists[0] = newlista
    
    
    newlist = np.array(lists[0])
            
    graphic = ['n']*cdna
    
    count = str(len(newlist))
    
    if int(count) == 0:
        print("Hmm.... There were no probes that fit the parameters. Check the full sequence field")
    else:
        df = pd.DataFrame()
        
        if BlastProbes == 'No':
            newlist = max33(maxprobe,newlist,numbr)
            newlist = newlist['reduced']
            count = str(len(newlist))
            print()
            print("From the given parameters, we were able to make "+count+" probe pairs.")
            print()
            print()
            a=0
            while a < len(newlist):
                seqs[a] = [newlist[a][0],str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]]),newlist[a][1]]
                graphic[newlist[a][0]:newlist[a][1]] = str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]])
                a+=1
            g = ''
            g = g.join(graphic)
            g = Seq(g)
            g = g.reverse_complement()

            df = pd.DataFrame.from_dict(seqs, orient='index', columns=['Start', 'seq', 'End'])
            df['InitiatorUp'] = upinit
            df['SpacerUp'] = uspc
            df['Probe1'] = df['seq'].str[27:52]
            df['Probe2'] = df['seq'].str[0:25]
            df['SpacerDw'] = dspc
            df['InitiatorDW'] = dninit
            df["Select"] = True
            columns = ["Select","InitiatorUp","SpacerUp","Start","Probe1","Probe2","End","SpacerDw","InitiatorDW"]
            df.drop(labels=['seq'], axis=1, inplace= True)
            df.index.name = "Pair"
            df = df[columns]

            grid = DataGrid(
                df,
                auto_fit_columns=True,
                auto_fit_params={'mode': 'all', 'max_rows': 2000, 'padding': 8},
                editable=True, 
                selection_mode="row",

                )
            grid.layout = {'width': '100%', 'height': '420px'}
            

            # # 1) Per-column renderers
            # renderers = {}

            # # (a) checkbox look for the boolean column
            # # We don't change the value; just how it's *drawn*
            # renderers['Select'] = TextRenderer(
            #     horizontal_alignment='center',
            #     # return '☑' for True, '☐' for False
            #     format=Expr("'☑' if cell.value else '☐'")
            # )

            # grid = DataGrid(df, editable=True, selection_mode="row", renderers= renderers)
            # grid.auto_fit_columns = True
            # display(grid)

            # print("Below is in IDT oPool submission_format.")
            # print("Copy and Paste the lines below into an XLSX file for submission to IDT starting from 'Pool name'.")
            # print()
            # print("Pool name, Sequence")

            # while a < len(newlist):
            #     seqs[a] = [newlist[a][0],str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]]),newlist[a][1]]
            #     graphic[newlist[a][0]:newlist[a][1]] = str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]])
            #     print(str(amplifier+'_'+name+'_'+count+'_Dla'+str(pause)+','+upinit+uspc+str(seqs[a][1][27:52])))
            #     print(str(amplifier+'_'+name+'_'+count+'_Dla'+str(pause)+','+str(seqs[a][1][0:25])+dspc+dninit))
            #     a+=1

        else:
            graphic = ['n']*cdna





        ## THE FOLLOWING SECTION CREATES A FASTA FILE FROM THE POTENTIAL PROBE SEQUENCES (BOTH 25BP PROBES COUPLED AS A SINGLE 52BP SEQUENCE INCLUDING A 2BP "nn" SPACER)        
            ## THE RESULTANT FASTA FILE IS BLASTED AGAINST THE USER SPECIFIED TRANSCRIPTOME FASTA 
            ## PROBES THAT MATCH A SEQUENCE IN BLAST WITH A LENGTH MATCH, 60BP > X > 40BP, AND AN E-VALUE < 1E-15 ARE KEPT, OTHERS ARE DISCARDED


        if BlastProbes == "Yes":
            print()
            print("BLASTn of probes in progress, this may take a few minutes.")
            print()
            seqs={} 
            remove = pd.DataFrame(columns = ["pos1","seq","pos2","fasta","num"])
            a=0
            tmpFA = open((str(name)+"PrelimProbes.fa"), "w")
            while a < len(newlist):
                nm = str('>'+str(a))
                seqs[a] = [newlist[a][0],str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]]),newlist[a][1],nm,a]
                graphic[newlist[a][0]:newlist[a][1]] = str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]])
                remove.loc[a,['pos1','seq','pos2','fasta','num']] = [newlist[a][0] , str(fullseq[newlist[a][0]:(newlist[a][0]+25)]+"nn"+fullseq[(newlist[a][0]+27):newlist[a][1]]), newlist[a][1],nm, a] 
                tmpFA.write(nm)
                tmpFA.write('\n')
                tmpFA.write(seqs[a][1])
                tmpFA.write('\n')
                a+=1
            tmpFA.close()
            g = ''
            g = g.join(graphic)
            g = Seq(g)
            g = g.reverse_complement()
            # same as seq w/ no blast
            seqs_sub = {k: v[:3] for k,v in seqs.items()}
            
            df = pd.DataFrame.from_dict(seqs_sub, orient='index', columns=['Start', 'seq', 'End'])
            df['InitiatorUp'] = upinit
            df['SpacerUp'] = uspc
            df['Probe1'] = df['seq'].str[27:52]
            df['Probe2'] = df['seq'].str[0:25]
            df['SpacerDw'] = dspc
            df['InitiatorDW'] = dninit
            df["Select"] = False
            columns = ["Select","InitiatorUp","SpacerUp","Start","Probe1","Probe2","End","SpacerDw","InitiatorDW"]
            df.drop(labels=['seq'], axis=1, inplace= True)
            df.index.name = "Pair"
            df = df[columns]
            
        ## Probe BLAST setup and execution from FASTA file prepared in previous step

            cline = bn(query = str(name)+"PrelimProbes.fa", subject = db, outfmt = 6, task = 'blastn-short') #this uses biopython's blastn formatting function and creates a commandline compatible command 
            stdout, stderr = cline() #cline() calls the string as a command and passes it to the command line, outputting the blast results to one variable and errors to the other

            ## From results of blast creating a numpy array (and Pandas database)
            dt = [(np.str_,8),(np.str_,40),(np.int32),(np.int32),(np.int32),(np.int32),(np.int32),(np.int32),(np.int32),(np.int32),(np.float64),(np.float64)]
            blastresult = (np.genfromtxt(io.StringIO(stdout),delimiter = '\t',dtype = dt))# "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore")

            blast_db = pd.DataFrame(blastresult) 

            ## Create a list of pairs that pass the criteria (good) and do not (bad)

            good = list(blast_db.loc[(blast_db['f11']>=75.0) & (blast_db['f10']<=float(1e-13))].f0)
            good = [int(n) for n in good] 
            # If there is more than two probes, it might be off target too right? mak it problematic for now.
            multi = {x for x in good if good.count(x) > 1}
            good = set(good)
            good -= multi      

            bad = list(blast_db.loc[(blast_db['f11']>=60.0) & (blast_db['f10']>float(1e-13))].f0)
            bad = [int(n) for n in bad]
            bad = set(bad)
            bad |= multi

            in_good = df.index.isin(good)
            in_bad  = df.index.isin(bad)


            # if there is not any good hits raise error
            if not good:
                print("Hmm.... There were no probes that fit the parameters specified within the FASTA file.   ")
                print()
                print()
                print("  Try increasing the length of homopolymers tolerated, or BLAST against a different FASTA file.")
                print()
                print("  If BLASTing a heterologous sequence, i.e. GFP, this error could be because the RNA doesn't exist in your species. ")
                print()
                print()
                raise ValueError("Empty list of probes that passes the criteria")

            # Add blast info to df
            df["Status"] = np.select(
                [in_good & in_bad, in_good & ~in_bad, ~in_good & in_bad],
                ["Problematic", "Good", "Bad"],
                default=np.nan  # or "Unknown"
            )

            # if there is more probe pairs than requested, filter them.
            if (int(count) > numbr):            
            # if there is enough good sequences, just use that filter and max33
                df_good = df.query("Status == 'Good'")
                if (len(df_good) > numbr):
                    # filter newlist with df.index the use that to filter df
                    newlist = newlist[df_good.index]
                    newlist = max33(maxprobe,newlist,numbr)
                    df = df.loc[newlist['keep_idx']]
                    newlist = newlist['reduced'] #optional?
                else:
                    diff = len(df_good) - numbr
                    df_bad = df.query("Status != 'Good'")
                    df_bad.sample(diff)
                    df = pd.concat([df_good,df_bad])

            grid = DataGrid(
                df,
                auto_fit_columns=True,
                auto_fit_params={'mode': 'all', 'max_rows': 2000, 'padding': 8},
                editable=True, 
                selection_mode="row",

                )
            grid.layout = {'width': '100%', 'height': '420px'}
                


        ## This loop takes the data from the blast result and filters out probe pairs that do not meet criteria
            ## by setting a length match requirement this eliminates off-target pairs and half-pairs
            ## the e-value threshold ensures that the probe is a good match to the target


            # print()
            # print()
            # print("The probes are provided below in IDT oPool submission format.")
            # print("Copy and Paste the lines below into an XLSX file for submission to IDT starting from 'Pool name'.")
            # print()
            # print()
            # print("Pool name, Sequence")
            # a=0
            # b=0
            # seqs1={}
            # while a < len(seqs):
            # #while a < len(uniques):
            #     tmp = (seqs[a])
            #     print(str(amplifier+'_'+name+'_'+count+'_Dla'+str(pause)+','+upinit+uspc+str(tmp[1][27:52])))
            #     print(str(amplifier+'_'+name+'_'+count+'_Dla'+str(pause)+','+str(tmp[1][0:25])+dspc+dninit))
            #     graphic[tmp[0]:tmp[2]] = str(tmp[1])
            #     seqs1[b]=tmp
            #     b+=1
            #     a+=1
            # seqs=seqs1
            
            
            g = ''
            g = g.join(graphic) 
            g = Seq(g)
            g = g.reverse_complement()



    
    
    
    

        return(grid)


    
    # output(cdna,g,fullseq,count,amplifier,name,pause,seqs)
    
    
        
    # print()
    # print()
    # if report == 'Yes':
    #     print("Run "+str(date.today())+"\n   with settings: \n\t5'Pause:\t"+str(pause)+" \n\tChoice of probe set:\t"+str((choose))+"\tPair used: "+str(choice)+" \n\tLength of acceptable polyA/polyT runs:\t"+str(polyAT)+" \n\tLength of acceptable polyC/polyG runs:\t"+str(polyCG)+" \n\tBLASTn of Probes:\t"+str((BlastProbes))+" \n\tRemoval of probes with low quality BLAST hits:\t"+str((dropout)) )

    # print()
    # print()
    # print("References: ")
    # print()
    # print(" See Choi et al. 2018 Development for HCR3.0 methodology details ")
    # print(" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6031405/ ")
    # print()
    # print(" See Wang et al. 2020 BioRxiv for expanded amplifier information ")
    # print(' https://www.biorxiv.org/content/10.1101/274456v3 ')