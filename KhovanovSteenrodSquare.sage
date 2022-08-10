"""
DOES NOT WORK WITH UNKNOT OR HOPF LINK COMPONENTS
"""


# Import packages needed to run the following code
from snappy import *
import numpy as np
import itertools
from collections import Counter

##################################################################################################
########################### Functions used in Link class construction  ###########################
##################################################################################################

def get_all0vertex(link):
    """
    This function defines the all 0s vertex of the n-cube diagram for its input, which is a Sage link object.
    It returns a list of cycles, defined using the PD code of this link.
    This function follows a similar process to get the complete smoothing at the all 0s vertex in Lipshitz and 
    Sarkar's KhovanovSteenrod.
    """
    pd = link.mirror_image().pd_code()
    circles=[]
    crossings=flatten(pd)
    travelled=[False,]*len(crossings)
    while False in travelled:
        temp=[]
        complete=False
        start=travelled.index(False)
        now = start
        while not complete:
            travelled[now]=True
            if now%2 == 0:
                temp = temp + [-1-int(now/4),]
                now=now+1;
            else:
                temp = temp + [1+int(now/4),]
                now=now-1;
            travelled[now]=True
            for j in range(len(crossings)):
                if crossings[j] == crossings[now]:
                    if j == start:
                        complete = True
                    if not travelled[j]:
                        now = j
        circles=circles+[temp,]
    return(circles)

def relabel(circles,flat, signed = False):
    max_val = max(flat)+1
    crossings = np.arange(max_val,max_val+len(set(flat))+1)
    new_circles = []
    for circle in circles:
        if signed == True:
            new_circle = [sign(val)*(crossings[abs(val)-1]) for val in circle]
        else:
            new_circle = [crossings[abs(val)-1] for val in circle]
        new_circles.append(new_circle)
    return(new_circles)

def generate_ncube(link, all_0_vertex):
    """
    This funtion takes in an oriented link along with its the all 0s vertex and defines the rest of this link's n-cube diagram.
    It returns a dictionary whose keys are a tuple of 0s and 1s and whose values are the corresponding complete smoothing of the link.
    This function follows a similar process to how the ncube or hypercube is generated in Lipshitz and Sarkar's KhovanovSteenrod.
    """
    ncube = {tuple([0 for ind in range(len(link.pd_code()))]): all_0_vertex}
    for num in range(1, len(link.pd_code())+1):
        for subset in itertools.combinations(np.arange(len(link.pd_code())),num):
            vertex = tuple([1 if ind in subset else 0 for ind in range(len(link.pd_code()))])
            last_vertex = tuple([1 if ind in subset[1:] else 0 for ind in range(len(link.pd_code()))])
            new_cycles = flip(ncube[last_vertex],last_vertex,subset[0]+1)[0]
            ncube[vertex] = new_cycles
    return(ncube)


##################################################################################################
###########################               The Link class               ###########################
##################################################################################################

class Link():
    """
    This class represents an oriented link.

    Inputs:
    - input_data: This is either the PD code or DT code for an oriented link. It is a list of lists.
    - input_type: The is either the string "PD_code" or "DT_code" depending on the value of input_data.
    - mirror_image: This is a Boolean, where its value is True if we should take the mirror image of
                    the inputted PD or DT code, and False otherwise. (default True)
    - print_progress: This is a Boolean that specifies whether or not we should print out progress
                      logs. (default True)

    Methods:
    - repr: This function returns a string representation describing this link.
    - _enhanced_states_comments: This function computes the enhanced states of this link.

    Data:
    - sage_link: the Sage link version of this link
    - pd_code: the PD code
    - writhe: the writhe of this link
    - orientation: a list of 1s and -1s describing whether each crossing is positive or negative
    - n: the number of crossings in this link diagram
    - n_plus: the number of positive crossings
    - n_minus: the number of negative crossings
    - enhanced_states: the enhanced states of this link
    - all_0_vertex: a description of the complete smoothing of this link that corresponds to the all 0s vertex of its n-cube diagram
    - ncube: the n-cube diagram for this link, used for computing its Khovanov homology
    - correspondences: the correspondences associated to this link via its n-cube diagram (F^q : 2^n -> B, the Burnside category.
                       See "KHOVANOV HOMOTOPY TYPE, BURNSIDE CATEGORY, AND PRODUCTS" by Lawson, Lipshitz, and Sarkar). This will be
                       initialized as an empty dictionary and filled as correspondences are computed for this link.
    """
    def __init__(self, input_data, input_type, mirror_image = False, print_progress = True):
        if print_progress:
            print("Initializing the link.")
        self.print_progress = print_progress
        self.input_data = input_data
        self.input_type = input_type
        self.mirror_image = mirror_image
        if input_type == "PD_code":
            self.sage_link = sage.knots.link.Link(input_data)
        elif input_type == "DT_code":
            self.sage_link = Knots().from_dowker_code(input_data)
        else:
            raise ValueError(f"Input type {input_type} is not supported. Please use either 'PD_code' or 'DT_code' as input type.")
        # For whatever reason, the code here has the PD code information, and I believe DT code information too, the wrong way around.
        # What we usually consider to be the PD code for a given link, this code considers as the PD code for the mirror image.
        # Thus, we apply a mirror image here by default to account for that.
        if not mirror_image:
            self.sage_link = self.sage_link.mirror_image()
        
        if print_progress:
            print("Computing the PD code.")
        self.pd_code = self.sage_link.pd_code()
        if print_progress:
            print("Computing the writhe.")
        self.writhe = self.sage_link.writhe()
        if print_progress:
            print("Computing the orientation.")
        self.orientation = self.sage_link.orientation()
        if print_progress:
            print("Computing the number of positive crossings, number of negative crossings, and total number of crossings.")
        self.n = len(self.orientation)
        self.n_minus = self.orientation.count(-1)
        self.n_plus = self.orientation.count(1)
        
        if print_progress:
            print("Computing the enhanced states.")
        self.enhanced_states = self.sage_link._enhanced_states()
        if print_progress:
            print("Defining the all 0s vertex.")
        self.all_0_vertex = get_all0vertex(self.sage_link)
        if print_progress:
            print("Defining the n-cube diagram.")
        self.ncube  = generate_ncube(self.sage_link, self.all_0_vertex)

        if print_progress:
            print("Initializing the dictionary that will hold the correspondences associated to this link.")
        self.correspondences = dict()


    def __repr__(self):
        """
        This function returns a string representation describing this link using its PD code.
        """
        return("Link with PD code: " + self.sage_link.pd_code())


    def _enhanced_states_comments(self):
        """
        --------------------------------------------------------------------------------------------------------------------
        DISCLAIMER: This function is not actually used in this script. It is a copy of the Sage function _enhanced_states(),
        which is used in our define_khovanov_complex() function (taken from Sage's _khovanov_homology_cached() function), but 
        with additional comments to help explain what the _enhanced_states() functions does.
        See the Sage function here: https://github.com/sagemath/sage/blob/f8df80820dc7321dc9b18c9644c3b8315999670b/src/sage/knots/link.py#L865
        --------------------------------------------------------------------------------------------------------------------
        
        This function generates all of the enhanced states (vertices and all possible labelings of the associated complete
        smoothing in the n-cube) for a given link.
        
        This function builds the complete smoothings associated to each vertex of the n-cube diagram for a link
        by iterating through each crossing in the original link and drawing a 0- or 1-smoothing, as indicated by
        the given vertex. To do this, we add one node for each strand in a crossing, where overstrands count as
        2 separate strands, and we add two nodes for each crossing, one corresponding to each pair of vertices
        we're about to join up. Then, we add 4 edges for each crossing, one connecting each strand node with its
        corresponding crossing node. Note that we label each crossing node by the triple containing the two strand
        nodes it connects to as well as the label for this crossing.
        
        Input:
            * self: the link we're interested in (as an object of class sage.knots.link.Link)
        Output:
            * tuple(states): all of the enhanced states of that link (as a tuple of tuples of the form:
                                (vertex: tuple of length n,
                                the cycles labeled by v_-: tuple of cycles, which are given by tuples of tuples,
                                the cycles labeled by v_+: tuple of cycles, which are given by tuples of tuples,
                                cohomological degree of this element: sage.rings.integer.Integer,
                                quantum degree of this element: sage.rings.integer.Integer)
                            )
        """
        
        writhe = self.writhe #number of positive crossings minus number of negative crossings
        crossings = self.pd_code #list of PD code labelings of each crossing
        ncross = len(crossings) #number of crossings
        smoothings = []
        nmax = max(flatten(crossings)) + 1 #used as first crossing label
        Gs = dict()
        print(f"There are {2**ncross} vertices to consider.")
        for i in range(2 ** ncross): #loop through all vertices of the n-cube
            v = Integer(i).bits()
            v = v + (ncross - len(v))*[0] #vertex of the n-cube
            G = Graph()
            for j, cr in enumerate(crossings): #loop through each coordinate of v and build the graph of the complete smoothing
                n = nmax + j #crossing label
                if not v[j]: #perform a 0-smoothing
                    # For negative crossings, we go from undercrossings to the left
                    G.add_edge((cr[3], cr[0], n), cr[0])
                    G.add_edge((cr[3], cr[0], n), cr[3])
                    G.add_edge((cr[1], cr[2], n), cr[2])
                    G.add_edge((cr[1], cr[2], n), cr[1])
                else: #perform a 1-smoothing
                    # positive crossings, from undercrossing to the right
                    G.add_edge((cr[0], cr[1], n), cr[0])
                    G.add_edge((cr[0], cr[1], n), cr[1])
                    G.add_edge((cr[2], cr[3], n), cr[2])
                    G.add_edge((cr[2], cr[3], n), cr[3])
            Gs[i] = [v,G]
            sm = set(tuple(sorted(x for x in b if isinstance(x, tuple)))
                    for b in G.connected_components(sort=False)) #label a cycle by the crossing nodes it goes through
            iindex = (writhe - ncross + 2 * sum(v)) // 2 #cohomological degree
            jmin = writhe + iindex - len(sm) #minimum quantum degree of a labeling of this vertex (doesn't get used)
            jmax = writhe + iindex + len(sm) #maximum quantum degree of a labeling of this vertex (doesn't get used)
            smoothings.append((tuple(v), sm, iindex, jmin, jmax)) #collect the cycles in the smoothing at v
        states = []  # we got all the smoothings, now find all the states
        for sm in smoothings: #loop through and collect all possible labelings of each smoothing & get its quantum degree
            for k in range(len(sm[1])+1):
                for circpos in itertools.combinations(sorted(sm[1]), k):  # Add each state
                    circneg = sm[1].difference(circpos)
                    #quantum degree (j) = writhe + cohomological deg + (# positive cycles) - (# negative cycles)
                    j = writhe + sm[2] + len(circpos) - len(circneg)
                    states.append((sm[0], tuple(sorted(circneg)), tuple(circpos), sm[2], j))
        return(tuple(states),Gs)




##################################################################################################
###############   Helper Functions from Lipshitz and Sarkar's KhovanovSteenrod    ################
###########    https://github.com/sucharit/KhovanovSteenrod/blob/master/main.sage     ############
##################################################################################################

def rotate(self, n):
    "Return result of rotating leftwards by n units."
    return(self[n:]+self[0:n])

def reflect(self):
    """ Returns the cycle (list) gotten by reflecting self and negating all the entries """
    n = len(self)
    return([-self[n-i-1] for i in range(len(self))])

def merge(self,other,n):
    "Merges cyclic lists self and other along the edge labeled n. Returns: CycList"
    if not (n in self):
        new_a = reflect(self)
    else:
        new_a = self
    if not (n in other):
        new_b = reflect(other)
    else:
        new_b = other
    a_loc = list(new_a).index(n)
    b_loc = list(new_b).index(n)
    answer = new_a[a_loc+1:]+new_a[0:a_loc]+[-n,]+new_b[b_loc+1:]+new_b[0:b_loc]+[-n,]
    return(answer)

def split(self,n):
    "Splits cyclic list self along the edge labeled n. Returns (b,c) of CycList's"
    adata = list(self)
    newn = n
    if not (n in adata):
        newn = -n
    first = adata.index(newn)
    second = first + adata[first+1:].index(newn) + 1
    b = adata[first+1:second]+[-newn,]
    c = adata[second+1:]+adata[0:first]+[-newn,]
    return(b,c)


def get_info(self, vertex):
    """
    This function is not directly from Lipshitz and Sarkar. It collects the data of their ResolvedKnot object in KhovanovSteenrod:
    https://github.com/sucharit/KhovanovSteenrod/blob/master/main.sage#L303-L312
    """
    N = len(vertex)
    if vertex != []:
        vertex = tuple(vertex)
    if vertex == []:
        vertex = self.N*(0,)
    edge_incidences = dict([(i+1,list()) for i in range(N)])
    for x in self:
        for y in x:
            edge_incidences[abs(y)].append(self.index(x))
    C = len(self)
    return(N,edge_incidences,C)


def flip(self, vertex, edge):
    "Returns the resolved link obtained by switching the resolution at edge. edge is a number from 1 to n."
    N,edge_incidences,C = get_info(self,vertex)
    new_vertex = vertex[0:edge-1]+((vertex[edge-1]+1)%2, )+vertex[edge:]
    if edge_incidences[edge][0]==edge_incidences[edge][1]: #This is the case of a split.
        sc = edge_incidences[edge][0] #This circle is being split.
        other_components = self[0:sc]+self[sc+1:]
        new_vertex_incidences = tuple(other_components) + split(self[sc],edge)
    if edge_incidences[edge][0]!=edge_incidences[edge][1]: #This is the case of a merge.
        mc1 = min((edge_incidences[edge][0],edge_incidences[edge][1]))
        mc2 = max((edge_incidences[edge][0],edge_incidences[edge][1]))
        other_components = self[0:mc1]+self[mc1+1:mc2]+self[mc2+1:]
        new_vertex_incidences = tuple(other_components) + (merge(self[mc1],self[mc2],edge),)
    return(new_vertex_incidences, new_vertex)

def second_index(data, elt):
    "Re,turns the second place that elt occurs in data"
    return data[data.index(elt)+1:].index(elt)+data.index(elt)+1

def ladybug(start, mid, edges):
    """
    start is a cyclist (initial vertex), mid is a list of 4 cyclists (middle 2 vertices, each of which has 2 cycles),
    and edges is a list of two edges, whose endpoints are linked S^0-s on start (1st and 2nd coords where start and final
    vertices differ).
    the 4 cyclists in mid are the 4 circles obtained by splitting start along the 2 edges
    start is in a ladybug configuration, so we want to output the matching on mid
    output is a pair ((self-matching on 0,1,2,3), 0 or 1) where 0 if matching was the right one, 1 if it was the wrong one.
    """
    answer=dict()

    a = abs(edges[0]) #first coord where start and final verts differ
    b = abs(edges[1]) #second coord where start and final verts differ
    
    #Rotate / reflect so that start looks like (a,?,-b,?,a,?,-b,?)
    if a in start:
        tempns = start#extra caution
    else:
        tempns = reflect(start)
    ns = rotate(tempns,tempns.index(a)) #ns = "new start"

    #Adjust the mids so that a and b (which are positive) both occur in all mids, and a is the first entry.
    nmids = list()
    for x in mid:
        if a in x:
            nmids.append(rotate(x,x.index(a)))
        else:
            temp = reflect(x)
            nmids.append(rotate(temp,temp.index(a)))            

    #Case that one of the positively oriented intervals is nonempty.
    if ns[1] != -b or ns[second_index(ns,a)+1] != -b:
        #Rotate some more, so that in (a,?_1,-b,?,a,?,-b,?), ?_1 is non-empty.
        if ns[1] == -b:
            ns = rotate(ns,second_index(ns,a))
        c = ns[1]
        #Now pair up the elts of mid
        matched_pair = [x for x in nmids if ( (x[1]==c) or (x[len(x)-1]==-c) )]
        answer[nmids.index(matched_pair[0])]=nmids.index(matched_pair[1])
        answer[nmids.index(matched_pair[1])]=nmids.index(matched_pair[0])
        unmatched = [i for i in range(4) if not (i in answer.keys())]
        answer[unmatched[0]]=unmatched[1]
        answer[unmatched[1]]=unmatched[0]
        return (answer, 0)

    #Case that both positively oriented intervals are empty.    
    else:
        #Rotate so that in (a,-b,?_1,a,-b,?_2), ?_1 is non-empty.
        if second_index(ns,a) == 2:
            ns = rotate(ns,second_index(ns,a)) #i.e., by 2...
        c = ns[2]
        matched_pair = [x for x in nmids if (x[(x.index(b)+1)%(len(x))]==c) or (x[len(x)-2]==-c)]#extra precaution with len
        answer[nmids.index(matched_pair[0])]=nmids.index(matched_pair[1])
        answer[nmids.index(matched_pair[1])]=nmids.index(matched_pair[0])
        unmatched = [i for i in range(4) if not (i in answer.keys())]
        answer[unmatched[0]]=unmatched[1]
        answer[unmatched[1]]=unmatched[0]
        return (answer,1)

def equiv(self, other):
    """
    This function has been changed slightly from Lipshitz and Sarkar's function of the same name.
    In Lipshitz and Sarkar's code, this function takes as input an instance of their CycList class,
    referring to a cyclic list, which is how they define cycles in complete smoothings of links. 
    We do not use that class here, so this function has been updated to account for us using lists 
    to define cycles.
    """
    md = list(self)
    od = list(other)
    nmd = list(reflect(self))
    md.sort()
    od.sort()
    nmd.sort()
    if md != od and nmd != od:
        return False
    else:
        #The sorted things could be equal without cyclists being equal. viz. the (2,n) torus links
        md = list(self)
        if nmd == od:
            md = list(reflect(self))
        od = list(other)
        #Modulo cyclic permutations, our md is now either same as od or same as reverse of od.
        if od[(od.index(md[0])+1)%(len(od))] == md[1%(len(od))]:#Just checks two consecutive elements.
            return True
        else:
            return False
        #This test still fails if we have an unknot or a Hopf link, or if a CycList has 1 element.



##################################################################################################
###############       Functions used in defining correspondences for a link       ################
##################################################################################################

def reverse_corr(corr):
    for key in corr:
        corr[key] = [list(reversed(tup)) for tup in corr[key]]
    return(corr)

def relabel(circles,flat, signed = False):
    max_val = max(flat)+1
    crossings = np.arange(max_val,max_val+len(set(flat))+1)
    new_circles = []
    for circle in circles:
        if signed == True:
            new_circle = [sign(val)*(crossings[abs(val)-1]) for val in circle]
        else:
            new_circle = [crossings[abs(val)-1] for val in circle]
        new_circles.append(new_circle)
    return(new_circles)

def get_correct_ladybug(start, mid, edges):
    """
    The ladybug computation from Lipshitz and Sarkar's code will either return the ladybug matching or the wrong matching.
    It labels whether the output is the ladybug matching, so if it is not, we will swap the matching so that it will be.
    """
    matching = ladybug(start, mid, edges)
    # If the matching is correct, return it as is.
    if matching[1] == 0:
        return(matching[0])
    # If the matching is correct, fix it and return that.
    else:
        fixed = dict()
        for key in matching[0].keys():
            if matching[0][key]%2 == 0: #0 or 2
                fixed[key] = matching[0][key]+1 #replace 0 by 1 and replace 2 by 3
            else: #1 or 3
                fixed[key] = matching[0][key]-1 #replace 1 by 0 and 3 by 2
        return(fixed)

def group_tuples(corr, link):
    """
    Used in define_correspondence() function.
    The input corr is a list of lists of potentially equivalent elements of the correspondence we are trying to define, 
    where these groups are looking for elements whose tuples are exactly equal to another tuple except in one coordinate,
    where that coordinate is not the first or last one since those must be equal. We will look through these groupings to
    determine which elements are actually equivalent.
    """
    ncube = link.ncube
    flat = remove_duplicates(flatten(link.pd_code))
    grouped = []
    for ind in corr:
        for group in corr[ind]:
            # If a group of elements has 2 elements in it, then these must be equivalent.
            if len(group) == 2:
                grouped.append(group)
            # If a group of elements has 4 elements in it, it's either a ladybug matching, 
            # in which case all 4 elements are truly equivalent, or it's not, in which case 
            # we really have 2 separate pairs of equivalent elements.
            elif len(group) == 4:
                pairs = [[],[]]
                vert1 = []
                vert2 = []
                vert1_negpos = []
                vert2_negpos = []
                for tup_ind, tup in enumerate(group):
                    if tup_ind == 0:
                        vert1.append(tup)
                        vert1_negpos.append([[set(lst2) for lst2 in lst] for lst in sage_to_ls(tup[ind])])
                    else:
                        if tup[ind][0] == vert1[0][ind][0]:
                            vert1.append(tup)
                            vert1_negpos.append([[set(lst2) for lst2 in lst] for lst in sage_to_ls(tup[ind])])
                        else:
                            vert2.append(tup)
                            vert2_negpos.append([[set(lst2) for lst2 in lst] for lst in sage_to_ls(tup[ind])])
                if (vert2[0][ind][0] != vert2[1][ind][0]):
                    raise Exception("While computing the ladybug matching, we found that one of the correspondences was not computed correctly.")
                if ((len(vert1) != 2) & (len(vert2) != 2)):
                    raise Exception("While computing the ladybug matching, we found that one of the correspondences was not computed correctly.")
                cycles1 = ncube[vert1[0][ind][0]]
                cycles2 = ncube[vert2[0][ind][0]]
                mids = determine_mids(cycles1, cycles2)
                all_start = ncube[group[0][ind-1][0]]
                edges = [ind2+1 for ind2,val in enumerate(group[0][ind-1][0]) if val != group[0][ind+1][0][ind2]]
                edges.sort()
                start = determine_start(cycles1, cycles2, all_start)
                matching = get_correct_ladybug(start,mids,edges)
                relabeled_mids = relabel(mids,flat)
                #now match up the enhanced states
                if set(relabeled_mids[0]) in vert1_negpos[0][0]: #first cycle in first vertex in vert1 is negative
                    match = matching[0]
                    for vert_ind, vert in enumerate(vert2):
                        if set(relabeled_mids[match]) in vert2_negpos[vert_ind][0]: #also negative
                            pairs[0] = [vert1[0],vert]
                            pairs[1] = [vert1[1],vert2[(vert_ind+1)%2]]
                elif set(relabeled_mids[0]) in vert1_negpos[0][1]: #first cycle in first vertex in vert1 is positive
                    match = matching[0]
                    for vert_ind, vert in enumerate(vert2):
                        if set(relabeled_mids[match]) in vert2_negpos[vert_ind][1]: #also positive
                            pairs[0] = [vert1[0],vert]
                            pairs[1] = [vert1[1],vert2[(vert_ind+1)%2]]
                else:
                    raise Exception("While computing the ladybug matching, we found that one of the correspondences was not computed correctly.")
                grouped.extend(pairs)
            # Under this first round of finding equivalent elements under a single face map delta^s_u, the groups should only either 
            # have 2 or 4 elements in them. Otherwise, throw and error.
            else:
                raise Exception("While computing the ladybug matching, we found that one of the correspondences was not computed correctly.")
    return(grouped)


def determine_mids(cycles1, cycles2):
    """
    This function is used in determining whether we have a ladybug matching.
    """
    cycles1_toadd = [cycle for cycle in cycles1 if cycle not in cycles2]
    cycles2_toadd = [cycle for cycle in cycles2 if cycle not in cycles1]
    if (len(cycles1_toadd) == 2) & (len(cycles1_toadd) == 2):
        return(cycles1_toadd + cycles2_toadd)
    cycles1_toremove = []
    cycles2_toremove = []
    for c1 in cycles1_toadd:
        for c2 in cycles2_toadd:
            if equiv(c1, c2):
                cycles1_toremove.append(c1)
                cycles2_toremove.append(c2)
    for cyc in cycles1_toremove:
        cycles1_toadd.remove(cyc)
    for cyc in cycles2_toremove:
        cycles2_toadd.remove(cyc)
    if (len(cycles1_toadd) == 2) & (len(cycles1_toadd) == 2):
        return(cycles1_toadd + cycles2_toadd)
    else:
        raise Exception

def determine_start(cycles1, cycles2, all_start):
    """
    This function is used in determining whether we have a ladybug matching.
    """
    cycles12 = cycles1 + cycles2
    start = list(copy(all_start))
    for c1 in all_start:
        for c2 in cycles12:
            if equiv(c1, c2) & (c1 in start):
                start.remove(c1)
    if len(start) == 1:
        return(start[0])
    print("we never found the starting cycle.....")

def sage_to_ls(es):
    """
    Sage stores enhanced states differently than Lipshitz and Sarkar do. In order to use Lipshitz and Sarkar's code for
    defining the ladybug matching, we transform the enhanced state from the Sage format to that  of Lipshitz and Sarkar.
    """
    neg = []
    for cycle in es[1]:
        c2 = [tup[2] for tup in cycle]
        neg.append(c2)
    pos = []
    for cycle in es[2]:
        c2 = [tup[2] for tup in cycle]
        pos.append(c2)
    return(neg,pos)

def join_groups(grouped):
    """
    Used in define_correspondence().
    Now that we have groups of equivalent elements under a single face map in the composition that is delta^s_U,
    we want to determine for the entire correspondence delta^s_U, which of these groups are equivalent to each other.
    """
    # First, collect together all groups of equivalent elements under a single face map that have the same first and last coordinate, 
    # meaning that they could possibly be equivalent under delta^s_U.
    possible_groupings = dict()
    elems_list = []
    for group in grouped:
        tup = tuple([group[0][0],group[0][-1]])
        if tup in elems_list: #check if an equivalent version of this element has already been included,
                                #and if so, add this to the list of equivalent versions of this element.
            possible_groupings[elems_list.index(tup)].append(group)
        else:
            elems_list.append(tup)
            possible_groupings[len(elems_list)-1] = [group]

    # For each of those groups of potentially equivalent elements, determine which are actually equivalent.
    # A pair of groups are equivalent if they have an element in common.
    corr = []
    for key in possible_groupings:
        accounted_for = []
        for pair1 in possible_groupings[key]:
            if pair1 not in accounted_for:
                accounted_for.append(pair1)
                to_add = True
                while to_add == True:
                    to_add = False
                    for pair2 in possible_groupings[key]:
                        if pair2 not in accounted_for:
                            # This line checks if the two groups we're considering as possibly equivalent have an element in common.
                            if any([x in pair1 for x in pair2]):
                                pair1.extend(pair2)
                                accounted_for.append(pair2)
                                to_add = True
                corr.append(pair1)
    return(corr)



##################################################################################################
###############                     The KhovanovComplex class                     ################
##################################################################################################

class KhovanovComplex():
    """
    This class represents the Khovanov complex for an oriented link. Its inputs are the same as for the Link class
    since the Khovanov complex must be for a paritcular link. We do not compute the cochain groups or differentials
    until a quantum grading is specified. Once we've done that, we can compute the homology of this complex, which
    is the Khovanov homology of the given oriented link.

    Inputs:
    - input_data: This is either the PD code or DT code for an oriented link. It is a list of lists.
    - input_type: The is either the string "PD_code" or "DT_code" depending on the value of input_data.
    - mirror_image: This is a Boolean, where its value is True if we should take the mirror image of
                    the inputted PD or DT code, and False otherwise. (default False)
    - print_progress: This is a Boolean that specifies whether or not we should print out progress
                      logs. (default True)

    Methods:
    - define_khovanov_complex: This function defines the Khovanov cochain complex and the associated differential matrices for a 
        given quantum degree.
    - define_pq_khovanov_complex: Return the cochain group C^(p,q) and the differentials with domain C^(p,q).
    - define_correspondence: This function defines the correspondence delta^s_U, where s is the augmented semi-simplicial degree,
        which in terms of the cohomological degree p is s = n - p - n_- - 1 (n = #crossings, n_- = #negative crossings).
    - get_st_positive_pairs: Find all (s,t)-positive pairs for all possible pairs (s,t) in delta^s_{U^-} x delta^s_{U^+}.
    - nabla_choose: This function helps to define the nabla_r map before we've applied the Z/2-linearization.
    - nabla_on_kh: Given the maximal, positive chains from Moran's sq^i construction, compute sq^i : H^(p,q) -> H^(p+i,q).
    - check_khovanov_nonzero: For the Steenrod square sq^i : H^(p,q) -> H^(p+i,q), check if both of these Khovanov homology groups 
        are nontrivial.
    - compute_sqi: This function computes sq^i: H^{p,q}(link;Z/2) -> H^{p+i,q}(link;Z/2).
    - get_degrees:  It determines which Steenrod operations might possibly be nonzero. To determine which Steenrod operations might 
        be nonzero, we compute the given link's Khovanov homology groups, and then find all pairs of nonzero groups in the same 
        quantum degree.
    - compute_all_sqi: This functions computes all Steenrod squares for a given link.
    - projection: This function projects an element from one correspondence to another.


    Data:
    - link: an instance of the Link class
    - cochain_groups: This is initialized as an empty dictionary and will hold the cochain groups in the Khovanov complex for the given link.
            The keys are quantum degrees q and the values are further dictionaries whose keys are pairs (p,q), where p is the cohomological degree
            and q is that same quantum degree.
    - differentials: This is initialized as an empty dictionary and will hold the differentials in the Khovanov complex for the given link.
            The keys are quantum degrees q and the values are further dictionaries whose keys are cohomological degrees p.
    - projection_storage: This is initialized as an empty dictionary and will hold projections between different correspondences for the given link.
            The keys are pairs of tuples, where the first tuple describes the domain correspondence in a given projection and the second tuple describes
            range correspondence. Each of these tuples has the form (p,q,U,order), just like the input to define_correspondence(), where p is the cohomological degree,
            q is the quantum degree, U is list describing the face map we're applying, and order is a Boolean specifying whether or not we care about how
            the order of the elements of U are applied as face maps. The value of this pair in the dictionary is a dictionary where the keys are indices 
            of elements in the domain correspondence and the value is the index of the element the key gets projected onto in the range correspondence.
    - max_chains: This is initialized as an empty dictionary and will hold the maximal, positive chains, as per Moran's sq^i construction, for the given link.
            The keys are triples (i,p,q), where i is the index of the sq^i that these maximal chains are part of the computation of, p is the cohomological
            degree, and q is the quantum degree. The values here are again dictionaries whose keys are pairs of indices (s,t) representing elements in
            delta_{U^-} x delta_{U^+}. The values here are pairs 
    - sqi: This is initialized as an empty dictionary and will hold the sq^i, as defined by Moran, on the Khovanov homology for the given link.
            The keys are triples (i,p,q), representing sq^i : H^(p,q)(L) -> H^(p+i,q)(L). The values represent this sq^i function. These values are also
            dictionaries whose keys are indices of the generators of H^(p,q)(L) and whose values are the image of that generator in H^(p+i,q)(L) as a vector.
    """
    def __init__(self, link, input_type, mirror_image = False, print_progress = True):
        if print_progress:
            print("Creating the link.")
        self.print_progress = print_progress
        if (input_type == "PD_code") | (input_type == "DT_code"):
            self.link = Link(link, input_type, mirror_image, print_progress)
        elif input_type == "link_class":
            self.link = link
        else:
            raise ValueError("The input_type must be either 'PD_code', 'DT_code', or 'link_class'.")
        
        if print_progress:
            print("Initializing the dictionary that will hold the cochain groups, differentials, projection maps, maximal chains, and sq^i functions.")
        self.cochain_groups = dict()
        self.differentials = dict()
        self.projection_storage = dict()
        self.max_chains = dict()
        self.sqi = dict()

    def define_khovanov_complex(self, q, ring=GF(2)):
        """
        --------------------------------------------------------------------------------------------------------------------
        DISCLAIMER: This function is almost an exact copy of the Sage function _khovanov_homology_cached(), except for the 
        exclusion of the last two lines of that function, which are:
            homologies = ChainComplex(complexes,check=False).homology()
            return tuple(sorted(homologies.items()))
        This Sage function computes the Khovanov homology of the given link at some particular quantum degree in all
        cohomological degrees by first defining the Khovanov cochain groups and differentials. Because we use the cochain
        elements and the differentials in our Steenrod square construction, we wanted to be able to return just these
        objects and not only the Khovanov homology groups. Thus, we removed the last two lines of the
        _khovanov_homology_cached() function from Sage, and kept the rest of the code the same. This is what is presented
        below, with added comments to help clarify what the code is doing and why.
        
        The other change I made was to replace the input parameter 'height' by 'q'. Since the code uses this parameter as the
        desired quantum degree that we are interested in and for us, height means something different from quantum degree, we
        use the variable 'q' here since that is what we typically use to refer to quantum degree.
        --------------------------------------------------------------------------------------------------------------------
        
        This function defines the Khovanov cochain complex and the associated differential matrices
        for a given quantum degree. This can be done for any coefficient ring, but we have set GF(2) as the default since
        that what we use for computing Steenrod operations.
        
        For the Khovanov cochain complex, we store just its generators as their enhanced states, gotten from the Sage function
        _enhanced_states(). Note that when we define the differential matrices, we first define the transpose of the matrix
        we want and then take the transpose of that matrix. I'm not sure why this is happening, but since this is how the Sage
        code was written, I left it as is.
        
        Inputs:
            * self: an instance of the KhovanovComplex class for the link we're interested in
            * q: quantum degree of the elements we're interested in (integer)
            * ring: coefficient ring (default: GF(2))
        
        Outputs: (also saved to the data of this KhovanovComplex instance)
            * cochain_groups: dictionary containing the generators of all of the Khovanov cochain groups in the given quantum degree,
                    where the keys are tuples (p,q) and the values are lists of all enhanced states in bigrading (p,q)
            * differentials: dictionary containing the differentials between cochain groups in the given quantum degree, where
                        the keys are the cohomological degree as an integer and the values are the differentials as Sage
                        matrices over the given ring, where the matrix could have 0 rows or 0 columns
        """
        # Check if these cochain groups and differentials have already been computed. If so, return them instead of recomputing them.
        if q in self.cochain_groups.keys():
            return(self.cochain_groups[q], self.differentials[q])
        
        crossings = self.link.pd_code #list of PD code crossing labels for all crossings
        ncross = len(crossings) #number of crossings
        states = [(_0, set(_1), set(_2), _3, _4) for (_0, _1, _2, _3, _4) in self.link.enhanced_states] #get all the enhanced states for this link
        cochain_groups = {}  #this will hold the enhanced states grouped by (p,q)
        for st in states:
            i, j = st[3], st[4] #i is cohomological degree, j is quantum degree
            if j == q: #store the enhanced state if its quantum degree matches what we're interested in
                if (i,j) in cochain_groups:
                    cochain_groups[i,j].append(st)
                else:
                    cochain_groups[i,j] = [st]
        differentials = {} #this will hold the differentials of the Khovanov cochain complex as matrices
        for (i, j) in cochain_groups:
            if (i+1, j) in cochain_groups:
                m = matrix(ring, len(cochain_groups[(i,j)]), len(cochain_groups[(i+1,j)])) #initialize the shape of the matrix with
                                                                        #len(cochain_groups[(i,j)]) rows and len(cochain_groups[(i+1,j)])) cols
                                                                        #NOTE: we'll transpose this matrix later
                for ii in range(m.nrows()):
                    V1 = cochain_groups[(i,j)][ii]
                    for jj in range(m.ncols()):
                        V2 = cochain_groups[(i+1, j)][jj]
                        V20 = V2[0] #get the vertex for this enhanced state
                        difs = [index for index,value in enumerate(V1[0]) if value != V20[index]] #coordinates where the two
                                                                                                #vertices differ
                        #the next line checks if (1) the vertices differ in exactly one coordinate and (2) the positive cycles
                        #of one enhanced state do not intersect the negative cycles of the other (in both possible ways), which
                        #checks that all cycles that are unchanged by this map do not change labelings.
                        if len(difs) == 1 and not (V2[2].intersection(V1[1]) or V2[1].intersection(V1[2])):
                            m[ii,jj] = (-1)**sum(V2[0][x] for x in range(difs[0]+1, ncross)) #define the differential pointwise
            else:
                m = matrix(ring, len(cochain_groups[(i,j)]), 0) #if the image space of this differential is empty, the matrix is empty
            differentials[i] = m.transpose() #transpose the differential matrix and store its value
            if not (i-1, j) in cochain_groups:
                differentials[i-1] = matrix(ring, len(cochain_groups[(i,j)]), 0) #store additional 0 matrices, presumably used for
                                                                    #computing Khovanov homoogy
        self.cochain_groups[q] = cochain_groups
        self.differentials[q] = differentials
        return(cochain_groups, differentials)

    def define_pq_khovanov_complex(self, p, q):
        """
        Return the cochain group C^(p,q) and the differentials with domain C^(p,q).
        """
        # Check if these cochain groups and differentials have already been computed. If so, return them instead of recomputing them.
        if q in self.cochain_groups.keys():
            if (p,q) in self.cochain_groups[q]:
                return(self.cochain_groups[q][(p,q)], self.differentials[q][p])
            # If q is in the keys of the cochain groups dictionary but (p,q) is not a key in that dictionary, 
            # then that means this cochain group was 0.
            else:
                return([], [])
        # If q is not in the keys of the cochain groups dictionary, then we haven't computed those cochain groups yet, so do so.
        else:
            gp, diff = self.define_khovanov_complex(q)
            return(self.define_pq_khovanov_complex(p,q))



    def define_correspondence(self, input_tuple):
        '''
        This function defines the correspondence delta^s_U, where s is the augmented semi-simplicial degree,
        which in terms of the cohomological degree p is s = n - p - n_- - 1 (n = #crossings, n_- = #negative crossings).
        An element of delta^s_U contains pairs of domain elements in X_s and range elements in X_{s-|U|}. However,
        here we represent these elements as the tuple of all enhanced states that the domain element passes through to
        get to the range element. If the order of U doesn't matter (order=False), we record all possible such tuples from
        the domain element to the range element for all orders of the elements of U. If the order of U does matter
        (order=True), we record only the tuple for each domain/range pair corresponding to applying |U|-many face maps
        starting at u_1, then u_2, and up to u_m. That is, we enforce that the order of the sequence U=[u_1, ..., u_m]
        matters, we return:
            delta^{s-(m-1)}_{u_m} circ delta^{s-(m-2)}_{u_{m-1}} circ ... circ
                                                                        delta^{s-1}_{u_2} circ delta^s_{u_1}.
                                                                        
        Note that for U=[u_1,u_2,...u_m], the elements of U tell us to change the u_i-th 0-coordinate from a 0 to 1, not
        that we should change the u_i-th overall coordinate from a 0 to 1, since for a given vertex, the u_i-th coordinate
        may already be a 1. This is why we use psi_inds to get the indices of the 0-coordinates that we should change to 1s,
        instead of just taking the u_i-th coordinate.
        
        This function has 3 parts:
            -Part 1: Collect all nontrivial elements of all correspondences from X_s to X_{s-|U|}.
            -Part 2: Keep only the elements in the correspondence delta^s_U. If order=True, in this part, we also
                    throw out any elements, as described above, gotten from applying the elements of U in a different order.
            -Part 3: Arbitrarily order the elements of this correspondence (used in determining (s,t)-goodness), and
                    if order=False, group equivalent elements from different orderings of U.
                    
        Inputs:
            * self: the link we're interested in
            * input_tuple is a tuple (p, q, U, order) where:
                * p: cohomological degree of the domain cochain group (integer)
                * q: quantum degree of the elements we're interested in (integer)
                * U: subscript on the face map we're applying (list)
                * order: whether or not we care about how the order of the elements of U are applied as face maps (Boolean)
        
        Output:
            * corr_set: collection of the elements of the correspondence delta^s_U, where the keys are natural numbers
                        providing an arbitrary order to the elements, and the values are either lists of equivalent elements
                        of the correspondence or a single element of the correspondence, depending on whether order is
                        True or False (dictionary)
        '''

        p = input_tuple[0]
        q = input_tuple[1]
        U = tuple(input_tuple[2])
        if len(input_tuple) > 3:
            order = input_tuple[3]
        else:
            order = False

        if input_tuple in self.link.correspondences.keys():
            return(self.link.correspondences[input_tuple])
        
        j=len(U) #tells us the number of cohomological degrees this correspondence increases us by
        if j == 0: #if U is the empty set, this correspondence is the identity
            corrs_set = {ind: [[elem,elem]] for ind,elem in enumerate(self.cochain_groups[q][(p,q)])}
            return(corrs_set) #return the identity if U is empty
        groups, diffs = self.define_khovanov_complex(q)
        if (p+j,q) not in groups:
            return(dict()) #return nothing if the range of this correspondence is empty
        if len([u for u, x in enumerate(groups[(p+j,q)][0][0]) if x == 1]) <= max(U):
            # Raise an exception if the correspondence is trying to change coordinates from 0 to 1 that don't exist
            # e.g., Suppose U = [5] and the vertices in C^(p,q) only have 3 0-coordinates. This correspondence
            # doesn't make sense since there is no 5-th 0-coordinate to change to a 1.
            raise Exception("You're trying to generate a nonsensical correspondence")
        U = list(U)
        U.reverse()
        
        #part 1: collect all elements in all delta^s_V, where |V| = |U|, namely V=[v_1, v_2, ..., v_m]
        for step in range(1,j+1): #define union of all face maps delta^s_V for all |V|=|U| via the composition:
                                #delta^{s-(m-1)}_{v_m} circ delta^{s-(m-2)}_{v_{m-1}} circ ... circ delta^s_{v_1}
            if step == 1: #compute delta^s_{v_1} for all v_1
                pairs1 = self.differentials[q][p].nonzero_positions() #each nonzero entry in position (i,j) in the differential matrix
                                                        #gives us an element of the correspondence:
                                                        #(j-th element of C^{p+1,q}, i-th element of C^{p,q})
                tuples = [[pair[::-1]] for pair in pairs1] #reverse the order of the pair as written in the line above to get:
                                                        #(i-th element of C^{p,q}, j-th element of C^{p+1,q})
            else: #post-compose delta^{s-(step-1)}_{v_step} with delta^s_{[v_1,v_2,...,v_{step-1}]}
                current_step = []
                nonzeros = diffs[p+step-1].nonzero_positions()
                nonzeros_rev = [pair[::-1] for pair in nonzeros] #as in the step 1 case, collect all pairs of elements of the
                                                                #correspondence delta^{s-(step-1)}_{v_step} and
                                                                #reverse their order
                for pair in itertools.product(tuples,nonzeros_rev): #compose the correspondences: match up pairs (a,b) & (b,c)
                    if (pair[0][-1][1] == pair[1][0]):
                        current_step.append([elem for elem in pair[0]]+[pair[1]])
                tuples = current_step    
        
        #part 2: keep only the elements in delta^s_U (either for all orderings of U or only ordered as given)
        corrs = []
        count = 0
        for pair in tuples: #loop through every correspondence element written as: (a,b), (b,c), (c,d),
                        #where a,b,c,d are the indices of the elements within their appropriate bases, and both join
                        #the tuples to get (a,b,c,d), and rewrite each coordinate as the actual basis element.
            if j == 1: #if U has a single element, then each correspondence element is already (a,b), so we don't have
                    #to combine the tuples, only get the actual basis elements.
                elem = [groups[(p+ind,q)][tup] for ind,tup in enumerate(pair[0])]
            else: #if U has more than one element, we have to combine the tuples and then get the actual basis element.
                elem = [groups[(p+ind,q)][tup[0]] for ind,tup in enumerate(pair)] + [groups[(p+j, q)][pair[-1][1]]]
            inds1 = [u for u, x in enumerate(elem[-1][0]) if x == 1] #get the indices of the 1-coordinates in the range vertex
            psi_inds = [inds1[u] for u in U] #subselect those indices from inds1 corresponding to U, e.g., the u_i-th 0-coord.
            if order == False: #keep all tuples in delta^s_U for all permutations of U, which means that the
                            #u_i-th 0-coordinates of the domain vertex becomes a 1 for all u_i in U (in any order
                            #across the tuple).
                if (sum([elem[0][0][u] for u in psi_inds])==0) & (sum([elem[-1][0][u] for u in psi_inds])==j):
                    corrs.append(elem)

                tracker = 0
                for label_ind, label in enumerate(elem): #loop each subelement of each element of the correspondence and keep
                                                        #only the elements that change the u_i-th 0-coordinates in order
                                                        #as given in U. we do this by counting the number of times that the
                                                        #(i+1)-st subelement changes the u_i-th 0-coord of the domain vertex,
                                                        #and only including the elements where this happens j+1 times.
                    if (sum([label[0][psi_inds[u]] for u in range(label_ind)])==label_ind) & \
                    (sum([label[0][psi_inds[u]] for u in range(label_ind,j)])==0):
                        tracker += 1
                if tracker == j+1: #we want j+1 and not j because each element of the correspondence is a (j+1)-tuple
                    count += 1
            else: #keep only those tuples corresponding to delta^s_U where first the u_1-st 0-coordinate becomes a 1,
                #then the u_2-nd 0-coord of the original vertex becomes a 1, etc.
                tracker = 0
                for label_ind, label in enumerate(elem): #loop each subelement of each element of the correspondence and keep
                                                        #only the elements that change the u_i-th 0-coordinates in order
                                                        #as given in U. we do this by counting the number of times that the
                                                        #(i+1)-st subelement changes the u_i-th 0-coord of the domain vertex,
                                                        #and only including the elements where this happens j+1 times.
                    if (sum([label[0][psi_inds[u]] for u in range(label_ind)])==label_ind) & \
                    (sum([label[0][psi_inds[u]] for u in range(label_ind,j)])==0):
                        tracker += 1
                if tracker == j+1: #we want j+1 and not j because each element of the correspondence is a (j+1)-tuple
                    corrs.append(elem)
        
        #part 3: turn the list of correspondence elements into a dictionary whose keys are integers, giving an arbitrary
        #but fixed order to the elements of this correspondence. if we don't care about the order of U, this will also
        #group together all ways to get from each domain vertex to its range vertex by applying face maps using different
        #permutations of U.
        corrs_set = dict()
        # if we care about the order of the elements of U in delta^s_U or if U has exactly one element in it, then we've only
        # saved exactly the elements we need and we do not need to group equivalent elements since there aren't any.
        if (order == True) | (len(U) == 1):
            for ind,corr in enumerate(corrs):
                corrs_set[ind] = [corr]
        # if we do care about the order of the elements of U in delta^s_U and U has more than one element in it, then we have
        # to determine which elements are equivalent.
        else:
            if len(corrs) == 0:
                return(dict())                    
            corrs2 = {ind:[] for ind in range(1, len(corrs[0])-1)}
            # loop over all of the middle coordinates in a tuple for a given element, find all other elements that are exactly the same
            # as our element except at the given coordinate. This should collect groups of elements that are potentially equivalent under
            # a single face map. Later we will separate any of these groupings where the elements are not all equivalent. We will also
            # join together other groups of elements that are equivalent but under a different face map.
            for ind in range(1, len(corrs[0])-1):
                elems_list = []
                for corr in corrs:
                    tup = tuple(corr[ind2] for ind2 in range(len(corr)) if ind2 != ind) #get all vertices in this tuple
                    if tup in elems_list: #check if an equivalent version of this element has already been included,
                                        #and if so, add this to the list of equivalent versions of this element.
                        if (corrs2[ind][elems_list.index(tup)][0][0] != corr[0]) | (corrs2[ind][elems_list.index(tup)][0][-1] != corr[-1]):
                            # raise an exception if the two elements we're claiming could be equivalent do not have the same first and last coordinate
                            raise Exception('Uh oh! We ran into an error defining this correspondence.')
                        corrs2[ind][elems_list.index(tup)].append(corr)
                    else:
                        elems_list.append(tup)
                        corrs2[ind].append([corr])
            # given these lists of potentially equivalent elements, now determine which are actually equivalent
            grouped = group_tuples(corrs2,self.link)
            corr_list = join_groups(grouped)
            corrs_set = {ind: corr_list[ind] for ind in range(len(corr_list))}
        # Initially I had defined these correspondences going the other direction. When I realized that they should be defined
        # as they are now, I took the lazy approach to fixing this and simply reversed the tuple that I had computed using the 
        # process above.
        final_corr = reverse_corr(corrs_set)
        self.link.correspondences[input_tuple] = final_corr
        return(final_corr)

    def get_st_positive_pairs(self, i, p, q, U_triple, all_s_t, corr_odd_tuple, corr_even_tuple):
        """
        Find all (s,t)-positive pairs for all possible pairs (s,t) in delta^s_{U^-} x delta^s_{U^+}.
        This is part of Moran's sq^i construction. Note that since (s,t)-positive pairs are by
        definition (s,t)-good, and so we skip the initial step of finding (s,t)-good pairs and 
        directly find the (s,t)-positive pairs.
        """
        n_minus = self.link.n_minus
        odd_cap_even = set(U_triple[1]).intersection(set(U_triple[2])) #this is the set U^- cap U^+

        # find all partitions of U^- cap U^+, which are the possible (s,t)-positive pairs.
        pairs = []
        for size1 in range(len(odd_cap_even)+1):
            for set1 in itertools.combinations(odd_cap_even, size1):
                leftover = odd_cap_even.difference(set1)
                for size2 in range(1,len(leftover)+1):
                    pairs.extend([[set(set1), set(set2)] for set2 in itertools.combinations(leftover,size2)])
        st_pos = {s_t: [] for s_t in all_s_t}

        if self.print_progress:
            print(f"There are {len(all_s_t)} (s,t) pairs to consider.")
        count = 0
        # Loop through each potential (s,t)-positive pair to determine whether it is, in fact, (s,t)-positive.
        for s_t in all_s_t:
            count += 1
            if (self.print_progress) & (count % 10 == 0):
                print(f"We're currently considering element #{count}.")
            source1 = get_source_target(s_t[0], self.define_correspondence(corr_odd_tuple), self.define_pq_khovanov_complex(p+i,q)[0], source=True)
            source2 = get_source_target(s_t[1], self.define_correspondence(corr_even_tuple), self.define_pq_khovanov_complex(p+i,q)[0], source=True)
            # We only care about (s,t)-positive pairs that are in the image of Delta^*.
            if source1 == source2:
                # each pair in pairs is a tuple (W'', W^{circ}).
                # we want to know if a given pair is (s,t)-good/positive for the given (s,t).
                for pair in pairs:
                    # in order to properly project along both lambda maps,
                    # we need to first biject from delta_U to delta_W'' circ delta_W^(circ) circ delta_{remaining}.
                    pair0 = list(pair[0])
                    pair1_wo_0 = list(pair[1].difference(pair[0]))
                    ordered_U = list(odd_cap_even.difference(pair[0].union(pair[1]))) + pair1_wo_0 +pair0
                    deg = p + i - len(ordered_U)
                    ordered_U.reverse()
                    delta_Upm_tuple = (deg,q,tuple(ordered_U),True)
                    s_prime = self.projection(s_t[0],corr_odd_tuple,delta_Upm_tuple)
                    t_prime = self.projection(s_t[1],corr_even_tuple,delta_Upm_tuple)
                    inds1 = pair0
                    inds1.reverse()
                    deg1 = p+i-len(pair0)
                    lambda_w_corr_tuple = (deg1,q,tuple(inds1),False)
                    lambda_w_s = self.projection(s_prime,delta_Upm_tuple,lambda_w_corr_tuple)
                    lambda_w_t = self.projection(t_prime,delta_Upm_tuple,lambda_w_corr_tuple)
                    # Check that lambda_W''(s) = lambda_W''(t), which is the first (s,t)-good requirement.
                    if (lambda_w_s == lambda_w_t):
                        inds2 = psi(pair0, pair1_wo_0)
                        inds2.reverse()
                        deg2 = p + i - len(pair[0].union(pair[1]))
                        lambda_w_w_corr_tuple = (deg2,q,tuple(inds2),False)
                        lambda_w_w_s = self.projection(s_prime,delta_Upm_tuple,lambda_w_w_corr_tuple)
                        lambda_w_w_t = self.projection(t_prime,delta_Upm_tuple,lambda_w_w_corr_tuple)
                        # the next line checks 2 things:
                        #     (1) that the projection is not null, and
                        #     (2) that the pair is (s,t)-positive (not negative), by comparing lambda_{W'',W^{circ}}(s) with lambda_{W'',W^{circ}}(t)
                        if (lambda_w_w_s != None) & (lambda_w_w_t != None) & \
                        ((((p+n_minus-1+(3*i))%2 == 0) & (lambda_w_w_s < lambda_w_w_t)) | \
                        (((p+n_minus-1+(3*i))%2 == 1) & (lambda_w_w_s > lambda_w_w_t))):
                            st_pos[s_t].append(pair)
        return(st_pos)   


    def nabla_choose(self, i, p, q):
        """
        This function helps to define the nabla_r map before we've applied the Z/2-linearization, constructed by Moran 
        in his paper "Higher Steenrod Squares for Khovanov Homology". That is, it tells us about cochain elements
        rather than cohomology elements.
        The result that gets returned from this function is a dictionary whose keys are pairs of indices for elements in 
        the target cochain group C^(p+i,q). The value of each dictionary is again a dictionary whose keys are indices of 
        elements in the source cochain group C^(p,q). The values here are collection of sets U in mathcal{P}_a(b) for 
        which the original key pair is an (s,t)-positive pair.
        """

        # To define this nabla_r map, we need to loop over the collection mathcal{P}_a(b), where b = p + i + n_minus,
        # where p is the cohomology degree we're interested in, i is the index of the sq^i we're defining, and 
        # n_minus is the number of negative crossings in this link. The numbers a and b come from the map nable^{a choose b}
        # that is used to define nabla_r. In particular, we're only interested in sets in mathcal{P}_a(b) that have 
        # exactly i even-indexed elements and exactly i odd-indexed elements. We call the group of even-indexed elements U^+
        # and the group of odd-indexed elements U^-.
        n_minus = self.link.n_minus
        # ints contains two copies of all numbers 0, ..., b. We do this because mathcal{P}_a(b) contains lists where 
        # integers can be repeated at most twice.
        ints = [int(np.floor(j/2)) for j in np.arange((p+i+n_minus)*2)]
        # Get all length 2i subsets of the list 'ints', where we removed duplicated subsets.
        all_U = remove_duplicates([pair for pair in itertools.combinations(ints,2*i)])
        # For each length 2i subset, split it into its even-indexed and odd-indexed elements.
        all_split_U = [split_U(U) for U in all_U]
        # Keep only subsets U in mathcal{P}_a(b) that have exactly i even-indexed elements and exactly i odd-indexed elements.
        U_Uminus_Uplus = [U_triple for U_triple in all_split_U if (len(U_triple[1]) == i) & (len(U_triple[2]) == i)]
        if self.print_progress:
            print(f'There are {len(U_Uminus_Uplus)} sequences U in mathcal(P)_{2*i}({p+i+n_minus}) that we will consider.')
        
        # Initialize dictionaries to hold all pairs (s,t) for which we will find the (s,t)-positive pairs and to hold the the maximal chain information.
        all_s_t = {U_triple[0]: [] for U_triple in U_Uminus_Uplus}
        max_chains = dict()
        for U_triple in U_Uminus_Uplus:
            if self.print_progress:
                print(f'We are currently considering {U_triple}.')

            # We want to define delta^r_{U^-} wedge delta^r_{U^+}.

            # Define delta^r_{U^-}, which is part of Moran's definition for nabla^{n choose q}
            corr_odd_tuple = (p,q,tuple(U_triple[1]),False)
            corr_odd = self.define_correspondence(corr_odd_tuple)
            # Define delta^r_{U^+}, which is part of Moran's definition for nabla^{n choose q}
            corr_even_tuple = (p,q,tuple(U_triple[2]),False)
            corr_even = self.define_correspondence(corr_even_tuple)
            
            # Collect all possible (s,t). Note that the first coordinates of s and t must be the same to account for 
            # Moran pre-composing this map with Delta_n to define the nabla_r map.
            all_s_t[U_triple[0]] = [prod for prod in itertools.product(corr_odd, corr_even) \
                                    if corr_odd[prod[0]][0][0] == corr_even[prod[1]][0][0]]
            
            # If U^- and U^+ have an empty intersection, then the wedge of their face maps is equal to
            # the product of their face maps, so let's compute that product.
            if len(set(U_triple[1]).intersection(set(U_triple[2]))) == 0:
                # Loop through each pair (s,t) in the product delta^r_{U^-} x delta^r_{U^+}.
                for elem in all_s_t[U_triple[0]]:
                    # The pair (s,t) defined in 'elem' is defined in terms of their index in our arbitrary ordering on corr_odd and corr_even.
                    # Instead, we want to record their index within the Khovanov cochain group that is either corr_odd's or corr_even's domain or range.
                    source1 = get_source_target(elem[0], corr_odd, self.define_pq_khovanov_complex(p+i,q)[0], source=True)
                    source2 = get_source_target(elem[1], corr_even, self.define_pq_khovanov_complex(p+i,q)[0], source=True)
                    target1 = get_source_target(elem[0], corr_odd, self.define_pq_khovanov_complex(p,q)[0], target=True)
                    target2 = get_source_target(elem[1], corr_even, self.define_pq_khovanov_complex(p,q)[0], target=True)
                    if source1 != source2:
                        raise Exception("Uh oh! s and t don't have the same source, but they should if they are to be in the image of Delta_n.")
                    # Store this pair (s,t) in the max_chain dictionary.
                    if (target1, target2) in max_chains.keys():
                        if source1 in max_chains[(target1, target2)]:
                            max_chains[(target1, target2)][source1].append([U_triple[0], [[set(), set()]]])
                        else:
                            max_chains[(target1, target2)][source1] = [[U_triple[0], [[set(), set()]]]]
                    else:
                        max_chains[(target1, target2)] = {source1: [[U_triple[0], [[set(), set()]]]]}
            
            # If U^- and U^+ have a non-empty intersection, then we'll actually have o compute their wedge product.
            else:
                # For each pair (s,t), determine the pairs (W'', W^(circ)), which are disjoint subsets of U_triple[0], that are (s,t)-positive
                st_pos = self.get_st_positive_pairs(i, p, q, U_triple, all_s_t[U_triple[0]], corr_odd_tuple, corr_even_tuple)

                # Initialize a dictionary to hold the positive, maximal chains.
                st_max = dict()
                # For each pair (s,t), determine the positive, maximal chains from its (s,t)-positive pairs computed above.
                for s_t in st_pos:
                    current_chains = []
                    for pair in st_pos[s_t]:
                        if len(pair[1]) == 1:
                            current_chains.append([pair])
                    if len(current_chains) > 0:
                        tracker = True
                        while tracker == True:
                            tracker, current_chains = get_chains(current_chains,st_pos[s_t], set(U_triple[1]).intersection(set(U_triple[2])))
                    else:
                        current_chains = []
                    if len(current_chains) > 0:
                        st_max[s_t] = current_chains
                
                for s_t in st_max:
                    for elem in st_max[s_t]:
                        # The pair (s,t) defined in 'elem' is defined in terms of their index in our arbitrary ordering on corr_odd and corr_even.
                        # Instead, we want to record their index within the Khovanov cochain group that is either corr_odd's or corr_even's domain or range.
                        source1 = get_source_target(s_t[0], corr_odd, self.define_pq_khovanov_complex(p+i,q)[0], source=True)                    
                        source2 = get_source_target(s_t[1], corr_even, self.define_pq_khovanov_complex(p+i,q)[0], source=True)
                        target1 = get_source_target(s_t[0], corr_odd, self.define_pq_khovanov_complex(p,q)[0], target=True)
                        target2 = get_source_target(s_t[1], corr_even, self.define_pq_khovanov_complex(p,q)[0], target=True)
                        if source1 != source1:
                            raise Exception("Uh oh! s and t don't have the same source, but they should if they are to be in the image of Delta_n.")
                        # Store this pair (s,t) in the max_chain dictionary.
                        if (target1, target2) in max_chains.keys():
                            if source1 in max_chains[(target1, target2)]:
                                max_chains[(target1, target2)][source1].append([U_triple[0], elem])
                            else:
                                max_chains[(target1, target2)][source1] = [[U_triple[0], elem]]
                        else:
                            max_chains[(target1, target2)] = {source1: [[U_triple[0], elem]]}

        # Store the positive, maximal chain information that we've just computed.
        self.max_chains[(i,p,q)] = max_chains
        return(max_chains)
   

    def nabla_on_kh(self, i, p, q):
        """
        Given the maximal, positive chains from Moran's sq^i construction, compute sq^i : H^(p,q) -> H^(p+i,q).
        """
        groups, differentials = self.define_khovanov_complex(q)
        max_chains = self.max_chains[(i,p,q)]
        sq = {}
        # If either differential needed for computing H^(p,q) is not in differentials, which would be the case if it were trivial, 
        # add in a trivial group of the right dimensions its place so that the homology computation below works.
        if p not in differentials.keys():
            if (p,q) in groups.keys():
                ncols = len(groups[(p,q)])
            else:
                ncols = 0
            if (p+1,q) in groups.keys():
                nrows = len(groups[(p+1,q)])
            else:
                nrows = 0
            differentials[p] = matrix(GF(2), nrows,ncols,0)
        if p-1 not in differentials.keys():
            if (p-1,q) in groups.keys():
                ncols = len(groups[(p-1,q)])
            else:
                ncols = 0
            if (p,q) in groups.keys():
                nrows = len(groups[(p,q)])
            else:
                nrows = 0
            differentials[p-1] = matrix(GF(2), nrows,ncols,0)

        hom = differentials[p].right_kernel() / differentials[p-1].column_space()
        hom_basis = hom.basis()
        if len(hom_basis) == 0:
            if self.print_progress:
                print("The domain homology groups are empty, so this Steenrod square is trivial.")
            self.sqi[(i,p,q)] = []
            return([])
        kh_h = [matrix(GF(2),len(hom.lift(hom_basis[ind])),1,
                    hom.lift(hom_basis[ind])) for ind in range(len(hom_basis))]


        # If either differential needed for computing H^(p+i,q) is not in differentials, which would be the case if it were trivial, 
        # add in a trivial group of the right dimensions its place so that the homology computation below works.
        if p+i not in differentials.keys():
            if (p+i,q) in groups.keys():
                ncols = len(groups[(p+i,q)])
            else:
                ncols = 0
            if (p+i+1,q) in groups.keys():
                nrows = len(groups[(p+i+1,q)])
            else:
                nrows = 0
            differentials[p+i] = matrix(GF(2), nrows,ncols,0)
        if p+i-1 not in differentials.keys():
            if (p+i-1,q) in groups.keys():
                ncols = len(groups[(p+i-1,q)])
            else:
                ncols = 0
            if (p+i,q) in groups.keys():
                nrows = len(groups[(p+i,q)])
            else:
                nrows = 0
            differentials[p+i-1] = matrix(GF(2), nrows,ncols,0)

        hom = differentials[p+i].right_kernel() / differentials[p+i-1].column_space()
        hom_basis = hom.basis()
        if len(hom_basis) == 0:
            if self.print_progress:
                print("The range homology groups are empty, so this Steenrod square is trivial.")
            self.sqi[(i,p,q)] = []
            return([])
        kh_hi = [matrix(GF(2),len(hom.lift(hom_basis[ind])),1,
                    hom.lift(hom_basis[ind])) for ind in range(len(hom_basis))]


        # Loop through each generator of the domain Khovanov homology group.
        # For each generator, collect the maximal chains whose target is that generator.
        # The maximal chains help us define Moran's nabla_r map, and so this section is 
        # equivalent to taking the dual of Moran's nabla_r map, thus giving us the value 
        # of sq^i on each generator. The image of each element is an element 
        # in the cochain group. In the next step, we'll pass to the quotient to define 
        # sq^i from cohomology group generators to cohomology group generators.
        for gen_ind, gen in enumerate(kh_h):
            sq[gen_ind] = [0 for i in range(len(groups[(p+i, q)]))]
            elems = [elem[0] for elem in gen]
            inds = [ind for ind,val in enumerate(elems) if val == 1]
            # each key in max_chains lists the source and target.
            # the value of each key in max_chains is the list of positive maximal chains.
            for key in max_chains:
                # checks if source of this element is in this generator tensored with itself
                if (key[0] in inds) & (key[1] in inds):
                    for target in max_chains[key]:
                        sq[gen_ind][target] += len(max_chains[key][target])
                        sq[gen_ind][target] = sq[gen_ind][target] %2
            if min(differentials[p+i].dimensions()) > 0:
                if max(differentials[p+i]*matrix(GF(2), len(sq[gen_ind]),1,sq[gen_ind]))[0] != 0:
                    raise Exception(f'Generator {gen_ind}: Uh oh! the output of the nabla map on generator {gen_ind} is not a cycle!')


        if len(sq) == 0:
            if self.print_progres:
                print("The Steenrod square is trivial.")
            self.sqi[(i,p,q)] = []
            return([])

        # Above, we defined a map that's almost sq^i: H^(p,q) -> H^(p+i,q), but actually, its range is C^(p+i,q).
        # To pass the image of this map to the quotient to get cohomology elements, we will want to solve the
        # equation A*X=B, where the columns of B are the image in C^(p+i,q) of sq^i, and A is a matrix whose 
        # first m columns are a basis for H^(p+i,q) and whose next m' columns are the differential d^(p+i-1).
        # When we solve for the matrix X, its top m rows will give us the linear combinations of cohomology classes
        # that are the image of sq^i.
        for ind,vec in enumerate(kh_hi): # set the first m columns of A to be a basis for H^(p+i,q)
            if ind == 0:
                A = vec
            else:
                A = A.augment(vec)
        if differentials[p+i-1].ncols() > 0:
            A = A.augment(differentials[p+i-1]) #set the next m' columns of A to be the differential whose image lands in C^(p+i,q)
        for ind,key in enumerate(sq): # set the columns of B to be the image of sq^i in C^(p+i,q)
            if ind == 0:
                B = matrix(GF(2), len(sq[key]),1,sq[key])
            else:
                B = B.augment(vector(GF(2),len(sq[key]),sq[key]))

        # To compute the Steenrod square, we solve for X in the expression A*X=B. If there is an error, 
        # we want to be able to print out a useful message. Thus, we use try/except here, so that in the 
        # case that something went wrong, the except ValueError section will print out the matrices that
        # we're looking at to help with debugging.
        try:
            X = A.solve_right(B) #solve for X in A*X=B
            #the first m rows of X tell us what cohomology classes make up the elements of B,
            #where m is the number of generators of kh^{h+i,q}
            relevant_X = X[:len(kh_hi)]
            sq_kh = {key: relevant_X.submatrix(row=0, col=ind, nrows=-1, ncols=1).transpose() for ind,key in enumerate(sq.keys())}
            for key in sq_kh:
                if key == 0:
                    sq_matrix = sq_kh[key].T
                else:
                    sq_matrix = sq_matrix.augment(sq_kh[key].T)
            sq_kh["rank"] = sq_matrix.rank()
            if self.print_progress:
                print("Done computing the Steenrod square")
            self.sqi[(i,p,q)] = sq_kh
            return(sq_kh)
        except ValueError:
            if self.print_progress:
                print("Something about this Steenrod square computation isn't working right. The result of the sq^i computation is not a collection of cochains, as it should be.")
            self.sqi[(i,p,q)] = {"A":A, "B":B, "rank":"unknown"}
            return({"A":A, "B":B, "rank":"unknown"})

    def check_khovanov_nonzero(self, i, p, q):
        """
        For the Steenrod square sq^i : H^(p,q) -> H^(p+i,q), check if both of these Khovanov homology groups are nontrivial.
        Return True if both groups are nontrivial, meaning that sq^i could potentially be nontrivial. Return False if at least
        one of these groups is trivial, meaning that sq^i must be the 0 map.

        Note: This function assumes that if any Khovanov homology groups have been stored in self.homology for a given q, 
        then all of the homology groups in that quantum degree have been computed and stored. That is, if (p,q) is not a key in
        self.homology[q], then we assume it is because the group H^(p,q) is trivial.
        """

        # First, consider the case where the Khovanov homology groups in quantum degree q have already been computed.
        domain_range = [False, False]
        if q in self.homology.keys():
            if p in self.homology[q].keys():
                if self.homology[q][p].dimension() > 0:
                    domain_range[0] = True
                else:
                    return(False)
            else:
                return(False)
            if p+i in self.homology[q].keys():
                if self.homology[q][p+i].dimension() > 0:
                    domain_range[1] = True
                else:
                    return(False)
            else:
                return(False)
        # For sum reason, the 'all' function wasn't working. The line below checks that both elements in domain_range are true.
        if sum(domain_range) == 2:
            return(True)
        
        # Now, consider the case where the Khovanov homology groups in quantum degree q have not been computed yet.
        groups, differentials = self.define_khovanov_complex(q)
        hom = differentials[p].right_kernel() / differentials[p-1].column_space()
        self.homology[q][p] = hom
        if hom.dimension() == 0:
            return(False)
        hom = differentials[p+i].right_kernel() / differentials[p+i-1].column_space()
        self.homology[q][p+i] = hom
        if hom.dimension() == 0:
            return(False)
        return(True)


    def compute_sqi(self, i, p, q):
        '''
        This function computes sq^i: H^{p,q}(link;Z/2) -> H^{p+i,q}(link;Z/2).
        '''
        if self.print_progress:
            print(f"Computing sq^{i}: H^{(p,q)} -> H^{(p+i,q)}.")
        if (i,p,q) in self.sqi.keys():
            if self.print_progress:
                print("We've already computed this Steenrod square.")
            return(self.sqi[(i,p,q)])

        if self.print_progress:
            print("Computing Khovanov cochain group and differentials")
        gp, diff = self.define_khovanov_complex(q)
        if self.print_progress:
            print("Computing the nabla choose map on cochain groups")
        max_chains = self.nabla_choose(i, p, q)
        if self.print_progress:
            print("Computing the Steenrod square")
        sqi = self.nabla_on_kh(i, p, q)
        if self.print_progress:
            print(f"The Steenrod square is {sqi}.")
        return(sqi)


    
    def get_degrees(self):
        '''
        This function is used in the compute_all_sq() function below. It determines which Steenrod operations might possibly
        be nonzero, and so we should compute them. To determine which Steenrod operations might be nonzero, we compute the
        given link's Khovanov homology groups, and then find all pairs of nonzero groups in the same quantum degree. For each
        pair, we record the domain Khovanov homology group's cohomological degree (p) and quantum degree (q), as well as the
        value i of the Steenrod operation Sq^i between these two groups, whichg is equal to the difference in the groups'
        cohomological degrees.
        
        Input:
            * link: the link we're interested in (as an object of class sage.knots.link.Link)
            
        Output:
            * ipq: collection of triples, one for each potentially nonzero Steenrod square, determined as described above,
                where each triple has the form (i, p, q), where i is the Steenrod square degree, p is the domain Khovanov
                homology group's cohomological degree, and q is the groups' quantum degree (list)
        '''
        ipq = []
        num_cycles = len(self.link.pd_code)+1
        n_plus = self.link.n_plus
        n_minus = self.link.n_minus
        deg_shift = n_plus - 2*n_minus
        min_q = -num_cycles + deg_shift
        max_q = num_cycles + deg_shift
        for q in range(min_q, max_q+1):
            if self.print_progress:
                print(f"Determining potentially nonzero (i,p,q) values for q={q}")
            groups, differentials = self.define_khovanov_complex(q)
            pqs = [(p,q) for p in list(groups.keys())]
            for pairs in itertools.combinations(pqs, 2):
                i = pairs[1][0] - pairs[0][0]
                if self.check_khovanov_nonzero(i, pairs[0][0], q) == True:
                    ipq.append((i,pairs[0][0], q))
        return(ipq)
                

    def compute_all_sqi(self):
        '''
        This functions computes all Steenrod squares for a given link. It does this by first finding all potentially nonzero
        Steenrod squares, using the function get_degrees(). Then for each of these potentially nonzero operations, we apply
        the compute_sqi() function to compute that Steenrod square and store that output.
        
        This function includes some print statements to show what it's working on since this function is often very slow.
        
        Input:
            * link: the link we're interested in (as an object of class sage.knots.link.Link)
            
        Output:
            * sqs: collection of all Steenrod squares for a given link, stored as a dictionary with keys being the triple
                (i,p,q), referring to Sq^i : H^{p,q}(link;Z/2) -> H^{p+i,q}(link;Z/2), and with values being the output of
                the compute_sqi() function, which are also dictionaries. These subdictionaries have keys being indices of the
                generators of H^{p,q} and values being vectors of length equal to the number of generators of H^{p+i,q},
                where the vector tells us the linear combination of the generators of H^{p+i,q} of the image under Sq^i
                of each generator in H^{p,q}.
        '''
        
        if self.print_progress:
            print('Determining which Steenrod squares to consider...')
        degs = self.get_degrees() #find the potentially nonzero Steenrod operations
        for deg in degs: #loop through each potentially nonzero Steenrod square
            if self.print_progress:
                print(f'Considering sq^{deg[0]}: H^{(deg[1],deg[2])} -> H^{(deg[1]+deg[0], deg[2])}...')
            sq = self.compute_sqi(deg[0],deg[1],deg[2]) #compute each Steenrod square
            if len(sq) > 0:
                if "rank" in sq.keys():
                    self.sqi[deg] = sq
        return(self.sqi)



    def projection(self, key, corr1_tuple, corr2_tuple):
        '''
        This function finds the image of 'key' under the projection from the correspondence 'corr1' to the
        correspondence 'corr2'. The correspondences 'corr1' and 'corr2' don't have to have the same domain space or
        range space, but 'corr1' will be a "longer" correspondence, meaning that corr1 = delta^s_U and
        corr2 = delta^t_V where |U| >= |V|. Recall that we represent these correspondences by the entire tuple of
        enhanced states from one element to its image along the correspodence. Thus, to compute the value of this
        projection on 'key', we figure out which part of the tuple of enhanced states that represents 'key' matches up
        with the tuples in 'corr2'. We do this by determining the indices of the tuples in 'key' with the same number of
        1s as the first and last elements in the tuples of 'corr2'. Then we look for the tuple in 'corr2' equal to that
        restricted version of 'key'. In essence, this projection is the result of restricting 'key' to only part of its
        entire tuple.
        
        Inputs:
            * key: index of the element in corr1 that we want to project (integer)
            * corr1: the domain correspondence of the projection we want to compute- output of define_correspondence()  (dictionary)
            * corr2: the range correspondence of the projection we want to compute- output of define_correspondence() (dictionary)
        
        Output:
            * key2: index of the image of 'key' under the projection from corr1 to corr2 (integer) - if 'key' has no image
                    under this projection or if corr2 is empty, return None
        '''
        corr1 = self.define_correspondence(corr1_tuple)
        corr2 = self.define_correspondence(corr2_tuple)


        if len(corr2) == 0:
            return(None) #if the image correspondence of our projection map is empty, return an empty projection of key
                
        #first, check if we've already computed this projection
        if (corr1_tuple,corr2_tuple) in self.projection_storage.keys():
            if key in self.projection_storage[(corr1_tuple,corr2_tuple)]:
                return(self.projection_storage[(corr1_tuple,corr2_tuple)][key])
        
        for ind,vertex in enumerate(corr1[0][0]): #find the part of key that we're projecting onto
            if sum(vertex[0]) == sum(corr2[0][0][0][0]): #find the vertex in key that has the same number of 1s as the
                                                        #first element in each tuple of elements of corr2.
                first_coord = ind
            if sum(vertex[0]) == sum(corr2[0][0][-1][0]): #find the vertex in key that has the same number of 1s as the
                                                            #last element in each tuple of elements of corr2.
                second_coord = ind
        if first_coord == second_coord:
            for elem,key2 in itertools.product(corr1[key],corr2):
                if corr2[key2][0][0] in elem:
                    self.projection_storage = add_to_projection_storage(key, key2, corr1_tuple, corr2_tuple, self.projection_storage)
                    return(key2)
        for key2 in corr2:
            for elem1, elem2 in itertools.product(corr1[key],corr2[key2]):
                if elem1[first_coord:second_coord+1] == elem2:
                    self.projection_storage = add_to_projection_storage(key, key2, corr1_tuple, corr2_tuple, self.projection_storage)
                    return(key2)
        self.projection_storage = add_to_projection_storage(key, None, corr1_tuple, corr2_tuple, self.projection_storage)
        raise Exception(f"we never found the matching element to project onto for {corr1[key]} onto elements that look like {corr2[0]}")



def add_to_projection_storage(key1, key2, corr1_tuple, corr2_tuple, projection_storage):
    """
    Function to store the projection of the element at index key1 in the correspondence defined by corr1_tuple
    to the element at index key2 in the correspondence defined by corr2_tuple.
    """
    if (corr1_tuple,corr2_tuple) in projection_storage.keys():
        projection_storage[(corr1_tuple,corr2_tuple)][key1] = key2
    else:
        projection_storage[(corr1_tuple,corr2_tuple)] = {key1: key2}
    return projection_storage



#helper function used in both get_st_positive_pairs() and nabla_choose()

def get_source_target(ind, corr, gp, source=False, target=False):
    '''
    This function is used in the functions get_st_good_pairs() and nabla_choose(). Given an element of a correspondence,
    this function returns the source or target of this element (depending on which of 'source' or 'target' is True), which
    is a cochain. Note that exactly one of source or target must be True. Having both be True will just return the source.
    
    Inputs:
        * ind: index of the element of the correspondence that we want to find the source or target of (integer)
        * corr: correspondence that the element we're interested in lives in (dictionary - output of define_correspondence())
        * gp: collection of generators of the Khovanov cochain group that's the image of either the source or target map
              (dictionary - output of define_khovanov_complex())
        * source: indicator of whether we want the source of the element ind (Boolean)
        * target: indicator of whether we want the target of the element ind (Boolean)
    
    Output:
        - if source == True:
            * gp.index(corr[ind][0][0]): index of the image in the source cochain group of our correspondence
                                         element ind (integer)
        if target == True:
            *gp.index(corr[ind][0][-1]): index of the image in the target cochain group of our correspondence
                                         element ind (integer)
    '''
    if source == True:
        return(gp.index(corr[ind][0][0]))
    if target == True:
        return(gp.index(corr[ind][0][-1]))






def psi(set1, set2):
    '''
    The map psi_V : {U in P_q(n) | U cap V = emptyset} -> P_q(n-p) is defined as follows:
                U = (u_1, u_2, ..., u_q) -> W = (w_1, w_2, ..., w_q); w_j = u_j - |{v in V | v < u_j}|
    Basically, we reindex U after removing V from the 'universe'.
    
    NOTE: We require that U and V (here, set1 and set2) are disjoint.
    
    In order to find (s,t)-good pairs, we need to define lambda^r_{W''} and lambda^r_{W'',W^{circ}} both of whose domains
    are delta^r_{U^- cap U^+}. The maps are both essentially projections onto different parts of this domain, but
    to make sure we correctly compute these projections, we need to reindex W'' and W^{circ}.
    
    The map lambda^r_{W''} projects onto the latter portion of delta^r_{U^- cap U^+} where we do the W'' portion of
    that face map. However, we'll have already done the rest of the face map, and so we have to reindex W'' to account for
    this (that is, since W'' tells us which 0-coordinates to change to 1s, and since the first part of this face map 
    already changed some of the 0-coordinates to 1s, we have to reindex W'' to make sure we change the correct
    coordinates). This psi map does this reindexing of W'' after removing the other indices in (U^- cap U^+)\ W''.
    
    Similarly, the map lambda^r_{W'', W^{circ}} projects onto the middle portion of delta^r_{U^- cap U^+} where we
    do the W^{circ} portion of that face map. This time, we'll have already done the part of delta^r_{U^- cap U^+}
    corresponding to all the indices outside of W'' and W^{circ}. Thus, again, we need to reindex W^{circ} after removing
    the other indices in (U^- cap U^+)\(W'' cup W^{circ}).
    
    Inputs:
        * set1: the subscript on the psi_V map, which contains the elements we want to "remove from the universe" (set)
        * set2: the input into the psi_V map, which contains the elements we want to reindex (set)
    
    Outputs:
        * reindexed: the elements of set2, after being reindexed after removing the elements set1 (list)
    '''
    
    reindexed = [w - len([u for u in set1 if u < w]) for w in set2] #compute psi_{set1}(set2)
    return(reindexed)

     

def remove_duplicates(lst):
    """
    Remove duplicate entries from a list.
    """
    return([ii for n,ii in enumerate(lst) if ii not in lst[:n]])

def split_Ubar(Ubar):
    """
    Split a set of unique elements into its even and odd-indexed elements, where index here is as defined by Moran.
    """
    #Ubar is the list of the set of elements of U (each element appears exactly once)
    odd = [u for ind,u in enumerate(Ubar) if (u+ind+1)%2==1]
    even = [u for ind,u in enumerate(Ubar) if (u+ind+1)%2==0]
    return([Ubar,set(odd),set(even)])

def split_U(U):
    """
    Given a set of integers U, where each integer may appear more than once, split the set into its even and odd-indexed elements,
    where index here is as defined by Moran.
    """
    # Ubar_split is the split into odd-indexed and even-indexed subsets of the unique elements of U.
    Ubar_split = split_Ubar(remove_duplicates(U))
    counts = Counter(U)
    for (key,value) in counts.items():
        if value > 1:
            # For repeated elements, add the other copy of this element to the partition that its copy is not in.
            if key in Ubar_split[1]:
                Ubar_split[2].add(key)
            else:
                Ubar_split[1].add(key)
    return([U,Ubar_split[1],Ubar_split[2]])


def get_chains(current_chains, possible_elems, U_pm):
    '''
    This function is used in the nabla_choose() function.
    For a given pair (s,t), it finds all, if any, maximal chains from the (s,t)-positive pairs (W'', W^{circ}), 
    which are disjoint subsets of (U^- cap U^+). Since part of a maximal chain can occur in multiple maximal chains
    it gets a little messy trying to make sure we collect all possible maximal chains.
    '''
    if len(current_chains[0][-1][1]) == len(U_pm):
        return(False, current_chains)
    next_chains = []
    # Consider the last element (A,B) in each chain (so far). For each element a in A, see if (A\a, B cup a) is in possible_elems,
    # meaning that it is a positive pair and should be added to the current chain.
    for chain in current_chains:
        chain[-1][0] = set(chain[-1][0])
        for a in chain[-1][0]:
            if [chain[-1][0].difference(set([a])), chain[-1][1].union(set([a]))] in possible_elems:
                chain_copy = chain.copy()
                chain_copy.append([chain[-1][0].difference(set([a])), chain[-1][1].union(set([a]))])
                next_chains.append(chain_copy)
        for a in U_pm.difference(chain[-1][0].union(chain[-1][1])):
            if [chain[-1][0], chain[-1][1].union(set([a]))] in possible_elems:
                chain_copy = chain.copy()
                chain_copy.append([chain[-1][0], chain[-1][1].union(set([a]))])
                next_chains.append(chain_copy)
    if len(next_chains) == 0:
        return(False, [])
    if len(next_chains[0][-1][1]) == len(U_pm):
        return(False, next_chains)
    return(True, next_chains)






###################################################################

def define_link(input_data, input_type, mirror_image = False, print_progress = True):
    globals()['Link'] = Link(input_data, input_type, mirror_image, print_progress)

def compute_khovanov_complex(input_data, input_type, mirror_image = False, print_progress = True):
    globals()['KhovanovComplex'] = KhovanovComplex(input_data, input_type, mirror_image, print_progress)

# compute_khovanov_complex([[5,6,2,1],[3,4,6,5],[1,2,4,3], [11,12,8,7],[9,10,12,11],[7,8,10,9]], "PD_code", True, False)
# tr2_sqi = KhovanovComplex.compute_sqi(2, -6, -14)
# tr2_all_sqi = KhovanovComplex.compute_all_sqi()