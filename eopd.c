/*
 * Main developer: Nico Van Cleemput
 * 
 * Copyright (C) 2014 Ghent University.
 * Licensed under the GNU GPL, read the file LICENSE.txt for details.
 */

/* This program reads plane triangulations from standard in and
 * checks whether each 4-tuple is contained in an extended outer planar disc.   
 * 
 * 
 * Compile with:
 *     
 *     cc -o eopd -O4 eopd.c
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>


#define MAXN 34            /* the maximum number of vertices */
#define MAXE (6*MAXN-12)    /* the maximum number of oriented edges */
#define MAXF (2*MAXN-4)      /* the maximum number of faces */
#define MAXVAL (MAXN-1)  /* the maximum degree of a vertex */
#define MAXCODELENGTH (MAXN+MAXE+3)
#define MAX_EOPD ((MAXF)*((MAXF)-1)*((MAXF)-2)/6) /*maybe too few: currently = #triples*/

#define INFI (MAXN + 1)

typedef int boolean;

#undef FALSE
#undef TRUE
#define FALSE 0
#define TRUE  1

typedef unsigned long long int bitset;

#define ZERO 0ULL
#define ONE 1ULL
#define EMPTY_SET 0ULL
#define SINGLETON(el) (ONE << (el))
#define IS_SINGLETON(s) ((s) && (!((s) & ((s)-1))))
#define HAS_MORE_THAN_ONE_ELEMENT(s) ((s) & ((s)-1))
#define IS_NOT_EMPTY(s) (s)
#define IS_EMPTY(s) (!(s))
#define CONTAINS(s, el) ((s) & SINGLETON(el))
#define CONTAINS_ALL(s, elements) (((s) & (elements)) == (elements))
#define ADD(s, el) ((s) |= SINGLETON(el))
#define ADD_ALL(s, elements) ((s) |= (elements))
#define UNION(s1, s2) ((s1) | (s2))
#define INTERSECTION(s1, s2) ((s1) & (s2))
//these will only work of the element is actually in the set
#define REMOVE(s, el) ((s) ^= SINGLETON(el))
#define REMOVE_ALL(s, elements) ((s) ^= (elements))
#define MINUS(s, el) ((s) ^ SINGLETON(el))
#define MINUS_ALL(s, elements) ((s) ^ (elements))
//the following macros perform an extra step, but will work even if the element is not in the set
#define SAFE_REMOVE(s, el) ADD(s, el); REMOVE(s, el)
#define SAFE_REMOVE_ALL(s, elements) ADD_ALL(s, elements); REMOVE_ALL(s, elements)


typedef struct e /* The data type used for edges */ {
    int start; /* vertex where the edge starts */
    int end; /* vertex where the edge ends */
    int rightface; /* face on the right side of the edge
                          note: only valid if make_dual() called */
    struct e *prev; /* previous edge in clockwise direction */
    struct e *next; /* next edge in clockwise direction */
    struct e *inverse; /* the edge that is inverse to this one */
    int mark, index; /* two ints for temporary use;
                          Only access mark via the MARK macros. */
    bitset vertices;

} EDGE;

EDGE *firstedge[MAXN]; /* pointer to arbitrary edge out of vertex i. */
int degree[MAXN];
bitset neighbourhood[MAXN];

EDGE *facestart[MAXF]; /* pointer to arbitrary edge of face i. */
int faceSize[MAXF]; /* pointer to arbitrary edge of face i. */
bitset faceSets[MAXF];

EDGE edges[MAXE];

static int markvalue = 30000;
#define RESETMARKS {int mki; if ((markvalue += 2) > 30000) \
       { markvalue = 2; for (mki=0;mki<MAXE;++mki) edges[mki].mark=0;}}
#define MARK(e) (e)->mark = markvalue
#define MARKLO(e) (e)->mark = markvalue
#define MARKHI(e) (e)->mark = markvalue+1
#define UNMARK(e) (e)->mark = markvalue-1
#define ISMARKED(e) ((e)->mark >= markvalue)
#define ISMARKEDLO(e) ((e)->mark == markvalue)
#define ISMARKEDHI(e) ((e)->mark > markvalue)

int nv;
int ne;
int nf;

bitset eopdFaces[MAX_EOPD]; //contains extension
bitset eopdVertices[MAX_EOPD]; //doesn't contain third vertex of extension
int eopdCount = 0;

//////////////////////////////////////////////////////////////////////////////

////////START DEBUGGING METHODS

void printFaces(){
    int i, j;
    for(i=0; i<nf; i++){
        fprintf(stderr, "%d) ", i+1);
        for(j=0; j<nv; j++){
            if(CONTAINS(faceSets[i], j)){
                fprintf(stderr, "%d ", j+1);
            }
        }
        fprintf(stderr, "\n");
    }
}

void printFaceTuple(bitset tuple){
    int i;
    fprintf(stderr, "Face tuple: ");
    for(i=0; i<nf; i++){
        if(CONTAINS(tuple, i)){
            fprintf(stderr, "%d ", i+1);
        }
    }
    fprintf(stderr, "\n");
}

void printFaceTupleFaces(bitset tuple){
    int i, j;
    for(i=0; i<nf; i++){
        if(CONTAINS(tuple, i)){
            fprintf(stderr, "%d) ", i+1);
            for(j=0; j<nv; j++){
                if(CONTAINS(faceSets[i], j)){
                    fprintf(stderr, "%d ", j+1);
                }
            }
            fprintf(stderr, "\n");
        }
    }
}

void printVertexTuple(bitset tuple){
    int i;
    fprintf(stderr, "Vertex tuple: ");
    for(i=0; i<nv; i++){
        if(CONTAINS(tuple, i)){
            fprintf(stderr, "%d ", i+1);
        }
    }
    fprintf(stderr, "\n");
}

////////END DEBUGGING METHODS

boolean greedyExtendEopdAndStore(bitset currentEopdVertices, bitset currentEopdFaces){
    int i = 0;
    
    while(i < ne &&
            !(CONTAINS_ALL(currentEopdVertices, edges[i].vertices) &&
              !CONTAINS(currentEopdFaces, edges[i].rightface) &&
              (INTERSECTION(currentEopdVertices, neighbourhood[edges[i].next->end]) == edges[i].vertices))){
        i++;
    }
    
    if(i==ne){
        //store the eOPD
        eopdFaces[eopdCount] = currentEopdFaces;
        eopdVertices[eopdCount] = currentEopdVertices;
        eopdCount++;
    } else {
        //extend
        greedyExtendEopdAndStore(
                UNION(currentEopdVertices, faceSets[edges[i].rightface]),
                UNION(currentEopdFaces, SINGLETON(edges[i].rightface)));
    }
}

boolean findEOPD_impl(bitset currentEopdVertices, bitset currentEopdFaces, bitset remainingFaces, EDGE *lastExtendedEdge){
    //first check whether this is a covering eOPD
    if(IS_NOT_EMPTY(INTERSECTION(currentEopdFaces, remainingFaces))){
        //store the eOPD
        greedyExtendEopdAndStore(currentEopdVertices, currentEopdFaces);
        return TRUE;
    }
    
    //otherwise try extending the eOPD
    EDGE *extension = lastExtendedEdge->next;
    
    if(INTERSECTION(currentEopdVertices, neighbourhood[extension->next->end]) ==
            extension->vertices){
            //face to the right of extension is addable
            if(findEOPD_impl(UNION(currentEopdVertices, faceSets[extension->rightface]),
                    UNION(currentEopdFaces, SINGLETON(extension->rightface)),
                    remainingFaces, extension)){
                return TRUE;
            }
    }
    
    extension = lastExtendedEdge->inverse->prev->inverse;
    
    if(INTERSECTION(currentEopdVertices, neighbourhood[extension->next->end]) ==
            extension->vertices){
            //face to the right of extension is addable
            if(findEOPD_impl(UNION(currentEopdVertices, faceSets[extension->rightface]),
                    UNION(currentEopdFaces, SINGLETON(extension->rightface)),
                    remainingFaces, extension)){
                return TRUE;
            }
    }
    
    return FALSE;
}

boolean findEOPD(bitset tuple){
    int i, j;
    //first we check the stored eOPD's
    for(i = 0; i < eopdCount; i++){
        bitset intersection = INTERSECTION(tuple, eopdFaces[i]);
        if(HAS_MORE_THAN_ONE_ELEMENT(intersection)){
            return TRUE;
        }
    }
    
    //then we try to find a new eOPD
    for(i = 0; i < nf; i++){
        if(CONTAINS(tuple, i)){
            //try to find a eOPD with face i as extension
            bitset remainingFaces = MINUS(tuple, i);
            //we use each edge once as a possible shared edge 
            EDGE *sharedEdge = facestart[i];
            for(j = 0; j < 3; j++){
                //construct initial eopd
                int neighbouringFace = sharedEdge->inverse->rightface;
                bitset currentEopdVertices = faceSets[neighbouringFace];
                bitset currentEopdFaces = UNION(SINGLETON(i), SINGLETON(neighbouringFace));
                if(findEOPD_impl(currentEopdVertices, currentEopdFaces, remainingFaces, sharedEdge->inverse)){
                    return TRUE;
                }
                sharedEdge = sharedEdge->next->inverse;
            }
        }
    }
    return FALSE;
}

boolean findUncoveredFaceTuple_impl(bitset tuple, bitset tupleVertices, int position, int size){
    if((position == nf) && (size < 4)){
        return FALSE;
    }
    if(size < 3){
        //just extend and continue
        if(findUncoveredFaceTuple_impl(tuple, tupleVertices, position+1, size)){
            return TRUE;
        }
        if(IS_EMPTY(INTERSECTION(tupleVertices, faceSets[position]))){
            if(findUncoveredFaceTuple_impl(UNION(tuple, SINGLETON(position)),
                    UNION(tupleVertices, faceSets[position]), position+1, size+1)){
                return TRUE;
            }
        }
    } else if(size == 3){
        //search for eOPD and if none found: go to 4-tuple
        if(findEOPD(tuple)){
            return FALSE;
        }
        //no eOPD found: extending tuple
        if(findUncoveredFaceTuple_impl(tuple, tupleVertices, position+1, size)){
            return TRUE;
        }
        if(IS_EMPTY(INTERSECTION(tupleVertices, faceSets[position]))){
            if(findUncoveredFaceTuple_impl(UNION(tuple, SINGLETON(position)),
                    UNION(tupleVertices, faceSets[position]), position+1, size+1)){
                return TRUE;
            }
        }
    } else {// size == 4
        //search for eOPD
        return !findEOPD(tuple);
    }
    //if we get here then all tuples extending the current tuple were covered
    return FALSE;
}

void constructInitialEopds(){
    int i;
    
    EDGE *sharedEdge = facestart[0];
    int neighbouringFace = sharedEdge->inverse->rightface;
    bitset currentEopdVertices = faceSets[neighbouringFace];
    bitset currentEopdFaces = UNION(SINGLETON(0), SINGLETON(neighbouringFace));
    greedyExtendEopdAndStore(currentEopdVertices, currentEopdFaces);
    bitset coveredFaces = eopdFaces[0];
    
    i = 1; //0 is covered
    while(i < nf && CONTAINS(coveredFaces, i)){
        i++;
    }
    if(i==nf){
        return; //all faces covered
    }
    sharedEdge = facestart[i];
    neighbouringFace = sharedEdge->inverse->rightface;
    currentEopdVertices = faceSets[neighbouringFace];
    currentEopdFaces = UNION(SINGLETON(i), SINGLETON(neighbouringFace));
    ADD_ALL(coveredFaces, eopdFaces[1]);
    
    i = nf - 1;
    while(i > 1 && CONTAINS(coveredFaces, i)){
        i--;
    }
    if(i==1){
        return; //all faces covered
    }
    sharedEdge = facestart[i];
    neighbouringFace = sharedEdge->inverse->rightface;
    currentEopdVertices = faceSets[neighbouringFace];
    currentEopdFaces = UNION(SINGLETON(i), SINGLETON(neighbouringFace));
}

boolean findUncoveredFaceTuple(){
    //reset counter for eOPD's
    eopdCount = 0;
    
    //start by constructing some eOPD's to exclude many tuples
    constructInitialEopds();
    
    return findUncoveredFaceTuple_impl(EMPTY_SET, EMPTY_SET, 0, 0);
}

//=============== Writing planarcode of graph ===========================

void writePlanarCodeChar(){
    int i;
    EDGE *e, *elast;
    
    //write the number of vertices
    fputc(nv, stdout);
    
    for(i=0; i<nv; i++){
        e = elast = firstedge[i];
        do {
            fputc(e->end + 1, stdout);
            e = e->next;
        } while (e != elast);
        fputc(0, stdout);
    }
}

void writePlanarCodeShort(){
    int i;
    EDGE *e, *elast;
    unsigned short temp;
    
    //write the number of vertices
    fputc(0, stdout);
    temp = nv;
    if (fwrite(&temp, sizeof (unsigned short), 1, stdout) != 1) {
        fprintf(stderr, "fwrite() failed -- exiting!\n");
        exit(-1);
    }
    
    for(i=0; i<nf; i++){
        e = elast = firstedge[i];
        do {
            temp = e->end + 1;
            if (fwrite(&temp, sizeof (unsigned short), 1, stdout) != 1) {
                fprintf(stderr, "fwrite() failed -- exiting!\n");
                exit(-1);
            }
            e = e->next;
        } while (e != elast);
        temp = 0;
        if (fwrite(&temp, sizeof (unsigned short), 1, stdout) != 1) {
            fprintf(stderr, "fwrite() failed -- exiting!\n");
            exit(-1);
        }
    }
}

void writePlanarCode(){
    static int first = TRUE;
    
    if(first){
        first = FALSE;
        
        fprintf(stdout, ">>planar_code<<");
    }
    
    if (nv + 1 <= 255) {
        writePlanarCodeChar();
    } else if (nv + 1 <= 65535) {
        writePlanarCodeShort();
    } else {
        fprintf(stderr, "Graphs of that size are currently not supported -- exiting!\n");
        exit(-1);
    }
    
}


//=============== Reading and decoding planarcode ===========================

EDGE *findEdge(int from, int to) {
    EDGE *e, *elast;

    e = elast = firstedge[from];
    do {
        if (e->end == to) {
            return e;
        }
        e = e->next;
    } while (e != elast);
    fprintf(stderr, "error while looking for edge from %d to %d.\n", from, to);
    exit(0);
}

/* Store in the rightface field of each edge the number of the face on
   the right hand side of that edge.  Faces are numbered 0,1,....  Also
   store in facestart[i] an example of an edge in the clockwise orientation
   of the face boundary, and the size of the face in facesize[i], for each i.
   Returns the number of faces. */
void makeDual() {
    register int i, sz;
    register EDGE *e, *ex, *ef, *efx;

    RESETMARKS;

    nf = 0;
    for (i = 0; i < nv; ++i) {

        e = ex = firstedge[i];
        do {
            if (!ISMARKEDLO(e)) {
                facestart[nf] = ef = efx = e;
                faceSets[nf] = EMPTY_SET;
                sz = 0;
                do {
                    ef->rightface = nf;
                    ADD(faceSets[nf], ef->end);
                    MARKLO(ef);
                    ef = ef->inverse->prev;
                    ++sz;
                } while (ef != efx);
                faceSize[nf] = sz;
                ++nf;
            }
            e = e->next;
        } while (e != ex);
    }
}

void decodePlanarCode(unsigned short* code) {
    /* complexity of method to determine inverse isn't that good, but will have to satisfy for now
     */
    int i, j, codePosition;
    int edgeCounter = 0;
    EDGE *inverse;

    nv = code[0];
    codePosition = 1;

    for (i = 0; i < nv; i++) {
        degree[i] = 0;
        neighbourhood[i] = SINGLETON(code[codePosition] - 1);
        firstedge[i] = edges + edgeCounter;
        edges[edgeCounter].start = i;
        edges[edgeCounter].end = code[codePosition] - 1;
        edges[edgeCounter].vertices = UNION(SINGLETON(i), SINGLETON(code[codePosition] - 1));
        edges[edgeCounter].next = edges + edgeCounter + 1;
        if (code[codePosition] - 1 < i) {
            inverse = findEdge(code[codePosition] - 1, i);
            edges[edgeCounter].inverse = inverse;
            inverse->inverse = edges + edgeCounter;
        } else {
            edges[edgeCounter].inverse = NULL;
        }
        edgeCounter++;
        codePosition++;
        for (j = 1; code[codePosition]; j++, codePosition++) {
            if (j == MAXVAL) {
                fprintf(stderr, "MAXVAL too small: %d\n", MAXVAL);
                exit(0);
            }
            ADD(neighbourhood[i], code[codePosition] - 1);
            edges[edgeCounter].start = i;
            edges[edgeCounter].end = code[codePosition] - 1;
            edges[edgeCounter].vertices = UNION(SINGLETON(i), SINGLETON(code[codePosition] - 1));
            edges[edgeCounter].prev = edges + edgeCounter - 1;
            edges[edgeCounter].next = edges + edgeCounter + 1;
            if (code[codePosition] - 1 < i) {
                inverse = findEdge(code[codePosition] - 1, i);
                edges[edgeCounter].inverse = inverse;
                inverse->inverse = edges + edgeCounter;
            } else {
                edges[edgeCounter].inverse = NULL;
            }
            edgeCounter++;
        }
        firstedge[i]->prev = edges + edgeCounter - 1;
        edges[edgeCounter - 1].next = firstedge[i];
        degree[i] = j;

        codePosition++; /* read the closing 0 */
    }

    ne = edgeCounter;

    makeDual();

    // nv - ne/2 + nf = 2
}

/**
 * 
 * @param code
 * @param length
 * @param file
 * @return returns 1 if a code was read and 0 otherwise. Exits in case of error.
 */
int readPlanarCode(unsigned short code[], int *length, FILE *file) {
    static int first = 1;
    unsigned char c;
    char testheader[20];
    int bufferSize, zeroCounter;
    
    int readCount;


    if (first) {
        first = 0;

        if (fread(&testheader, sizeof (unsigned char), 13, file) != 13) {
            fprintf(stderr, "can't read header ((1)file too small)-- exiting\n");
            exit(1);
        }
        testheader[13] = 0;
        if (strcmp(testheader, ">>planar_code") == 0) {

        } else {
            fprintf(stderr, "No planarcode header detected -- exiting!\n");
            exit(1);
        }
        //read reminder of header (either empty or le/be specification)
        if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
            return FALSE;
        }
        while (c!='<'){
            if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
                return FALSE;
            }
        }
        //read one more character
        if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
            return FALSE;
        }
    }

    /* possibly removing interior headers -- only done for planarcode */
    if (fread(&c, sizeof (unsigned char), 1, file) == 0) {
        //nothing left in file
        return (0);
    }

    if (c == '>') {
        // could be a header, or maybe just a 62 (which is also possible for unsigned char
        code[0] = c;
        bufferSize = 1;
        zeroCounter = 0;
        code[1] = (unsigned short) getc(file);
        if (code[1] == 0) zeroCounter++;
        code[2] = (unsigned short) getc(file);
        if (code[2] == 0) zeroCounter++;
        bufferSize = 3;
        // 3 characters were read and stored in buffer
        if ((code[1] == '>') && (code[2] == 'p')) /*we are sure that we're dealing with a header*/ {
            while ((c = getc(file)) != '<');
            /* read 2 more characters: */
            c = getc(file);
            if (c != '<') {
                fprintf(stderr, "Problems with header -- single '<'\n");
                exit(1);
            }
            if (!fread(&c, sizeof (unsigned char), 1, file)) {
                //nothing left in file
                return (0);
            }
            bufferSize = 1;
            zeroCounter = 0;
        }
    } else {
        //no header present
        bufferSize = 1;
        zeroCounter = 0;
    }

    if (c != 0) /* unsigned chars would be sufficient */ {
        code[0] = c;
        if (code[0] > MAXN) {
            fprintf(stderr, "Constant MAXN too small: %d > %d \n", code[0], MAXN);
            exit(1);
        }
        while (zeroCounter < code[0]) {
            code[bufferSize] = (unsigned short) getc(file);
            if (code[bufferSize] == 0) zeroCounter++;
            bufferSize++;
        }
    } else {
        readCount = fread(code, sizeof (unsigned short), 1, file);
        if(!readCount){
            fprintf(stderr, "Unexpected EOF.\n");
            exit(1);
        }
        if (code[0] > MAXN) {
            fprintf(stderr, "Constant MAXN too small: %d > %d \n", code[0], MAXN);
            exit(1);
        }
        bufferSize = 1;
        zeroCounter = 0;
        while (zeroCounter < code[0]) {
            readCount = fread(code + bufferSize, sizeof (unsigned short), 1, file);
            if(!readCount){
                fprintf(stderr, "Unexpected EOF.\n");
                exit(1);
            }
            if (code[bufferSize] == 0) zeroCounter++;
            bufferSize++;
        }
    }

    *length = bufferSize;
    return (1);
}

//====================== USAGE =======================

void help(char *name) {
    fprintf(stderr, "The program %s checks extended outer planar discs in plane triangulations.\n\n", name);
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "\nThis program can handle graphs up to %d vertices.\n\n", MAXN);
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
}

void usage(char *name) {
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
         {"help", no_argument, NULL, 'h'}
    };
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case '?':
                usage(name);
                return EXIT_FAILURE;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }

    unsigned long long int numberOfGraphs = 0;
    unsigned long long int numberOfUncoveredGraphs = 0;

    /*=========== read planar graphs ===========*/

    unsigned short code[MAXCODELENGTH];
    int length;
    while (readPlanarCode(code, &length, stdin)) {
        decodePlanarCode(code);
        if(findUncoveredFaceTuple()){
            writePlanarCode();
            numberOfUncoveredGraphs++;
        }
        numberOfGraphs++;
    }
    
    fprintf(stderr, "Read %llu graph%s.\n", numberOfGraphs, 
                numberOfGraphs==1 ? "" : "s");
    fprintf(stderr, "Written %llu uncovered graph%s.\n", numberOfUncoveredGraphs, 
                numberOfUncoveredGraphs==1 ? "" : "s");
}
