import numpy as np
class Local:
    def __init__(self,str1,str2):
        self.flipped = False
        if(len(str1)>len(str2)):
            self.str1 = str2
            self.str2 = str1
            self.flipped = True
        else:
            self.str1 = str1
            self.str2 = str2
        #creating matrices
        main_matrix = np.zeros((len(self.str1)+1,len(self.str2)+1))
        match_check_matrix = np.zeros((len(self.str1),len(self.str2)))
        #Scores for match,mismatch and gap
        self.match = 1
        self.penalty = 0
        self.gap_penalty = -1

        self.cols = len(str2)
        #fillling up the match)check_matrix
        for i in range(len(self.str1)):
            for j in range(len(self.str2)):
                if(self.str1[i]==self.str2[j]):
                    match_check_matrix[i][j] = self.match
                else:
                    match_check_matrix[i][j] = self.penalty
        self.match_check_matrix = match_check_matrix
        #filling up the main matrix using Needleman-Wunsch
        #step 1 : initialization
        for i in range(len(self.str1)+1):
            main_matrix[i][0] = i* self.gap_penalty
        for j in range(len(self.str2)+1):
            main_matrix[0][j] = j* self.gap_penalty

        #step2: matrix filling
        for i in range(1,len(self.str1)+1):
            for j in range(1,len(self.str2)+1):
                main_matrix[i][j] = max(main_matrix[i-1][j-1]+match_check_matrix[i-1][j-1],
                                        main_matrix[i-1][j]+self.gap_penalty,
                                        main_matrix[i][j-1]+self.gap_penalty)
        self.main_matrix = main_matrix

    def align(self):
        #step 3: Backtracking
        lst_alignments = []
        #lst_path = []
        aligned_1 = ""
        aligned_2 = ""
        ti = len(self.str1)
        tj = len(self.str2)
        while(ti>0 and tj>0):
            score_current = self.main_matrix[ti][tj]
            score_diagonal = self.main_matrix[ti-1][tj-1]
            score_up = self.main_matrix[ti-1][tj]
            score_left = self.main_matrix[ti][tj-1]

            if(score_current == score_diagonal+self.match_check_matrix[ti-1][tj-1]):
                aligned_1 = self.str1[ti-1] + aligned_1
                aligned_2 = self.str2[tj-1] + aligned_2
                ti = ti - 1
                tj = tj - 1
              #  lst_path.append((ti,tj))
            elif(score_current == score_up+ self.gap_penalty):
                aligned_1 = self.str1[ti-1] + aligned_1
                aligned_2 = "-" + aligned_2
                ti = ti -1
             #   lst_path.append((ti,tj))
            elif(score_current==score_left+self.gap_penalty):
                aligned_1 = "-" + aligned_1
                aligned_2 = self.str2[tj-1]+aligned_2
                tj = tj - 1
            #    lst_path.append((ti,tj))

        if(self.flipped == True):
            temp = aligned_1
            aligned_1 = aligned_2
            aligned_2 = temp
        lst_alignments.append([aligned_1,aligned_2])
        #print(lst_alignments)
        # print(lst_path)
        # print(self.main_matrix)

    def makeGraph(self):
        graph = {}
        for i in range(1, len(self.str1)+1):
            # graph[(i,0)] = [(i-1,0)]
            # graph[(0,i)] = [(0,i-1)]
            for j in range(1,len(self.str2)+1):
                graph[(i,j)] = []
                score_current = self.main_matrix[i][j]
                score_diagonal = self.main_matrix[i - 1][j - 1]
                score_up = self.main_matrix[i - 1][j]
                score_left = self.main_matrix[i][j - 1]
                if( score_current == score_diagonal+self.match_check_matrix[i-1][j-1]):
                    graph[(i,j)] += [(i-1,j-1)]
                if(score_current == score_left + self.gap_penalty):
                    graph[(i,j)] += [(i,j-1)]
                if(score_current == score_up + self.gap_penalty):
                    graph[(i, j)] += [(i - 1, j)]
        return graph

    def find_AllPaths(self,graph,start,end,path = []):
        paths = []
        path = path + [start]
        if start == end:
            return [path]
        if start not in graph:
            return []
        for node in graph[start]:
            if node not in path:
                newpaths = self.find_AllPaths(graph,node,end,path)
                for newpath in newpaths:
                    paths.append(newpath)
        return paths

    def alignments(self):
        graph = self.makeGraph()
        tracks = self.find_AllPaths(graph,(len(self.str1),len(self.str2)),(0,0))
        alignments = []
        for track in tracks:
            alignment1 = ""
            alignment2 = ""
            last_step = (len(self.str1),len(self.str2))
            for step in track:
                i,j = last_step

                # if i == step[0] and j == step[1]:
                #     alignment1 = "-"+alignment1
                #     alignment2 = self.str2[j-1]+alignment2
                #     alignment1 = self.str1[i - 1] + alignment1
                #     alignment2 = "-" + alignment2
                if i == step[0]:
                    alignment1 = "-"+alignment1
                    alignment2 = self.str2[j-1]+alignment2

                elif j == step[1]:
                    alignment1 = self.str1[i-1] + alignment1
                    alignment2 = "-"+alignment2
                else:
                    alignment1 = self.str1[i-1]+ alignment1
                    alignment2 = self.str2[j-1]+ alignment2

                last_step = step
            alignment1 = alignment1[:-1]
            alignment2 = alignment2[:-1]
            if self.flipped == True:
                alignments.append([alignment2,alignment1])
            else:
                alignments.append([alignment1,alignment2])

        print(alignments)




NW = Local("SCYTHE","SCTHE")
NW.alignments()

NW = Local("AVNCCEGQHI","ARNDEQ")
NW.alignments()

NW = Local("SEQWENCE","SEQWENCE")
NW.alignments()

NW = Local("CYVPST","WYVPST")
NW.alignments()

NW = Local("WRYVPST","WYVPSAT")
NW.alignments()

NW = Local("ADMINS","ADMIRES")
NW.alignments()



