/*
 *  File      : duffybamtools.c
 *
 *  Created on: 17.06.2011
 *  Author    : W. Kaisers
 *  Content   : Function definitions in C for R package rbamtools
 *
 */

#ifndef duffybamtools_c
#define duffybamtools_c
#include "duffybamtools.h"

#define	BASE_1			1
#define	BAMTAG_MISMATCH		"MD"
#define MAX_CIGAR_PER_ALIGN  	5
#define MAX_MISMATCH_PER_ALIGN  10


static R_INLINE void clear_buf(char *c,unsigned n)
{
	int i;
	for(i=0;i<n;++i)
		c[i]=(char)0;
}

static R_INLINE void set_flag(bam1_t *align,_Bool val,unsigned pattern)
{
	if(val)
		align->core.flag=(align->core.flag) | pattern;
	else
		align->core.flag=(align->core.flag) & !pattern;
}

static R_INLINE int cigar2str(char *c,const bam1_t *align)
{
	if(align==NULL) 
		return 0;

	uint32_t len=align->core.n_cigar;
	if(len==0) {
		strcpy( c, "*");
		return 1;
	}
	uint32_t *cigar=bam1_cigar(align);
	char buf[128];

	sprintf(buf,"%lu",(unsigned long) (cigar[0] >> BAM_CIGAR_SHIFT));
	strcpy(c,buf);
	if((cigar[0]&BAM_CIGAR_MASK)>=strlen(CIGAR_TYPES))	// Error
		return 0;
	strncat(c,&(CIGAR_TYPES[cigar[0] & BAM_CIGAR_MASK]),1);


	uint32_t i;
	for(i=1;i<len;++i)
	{
		sprintf(buf,"%lu",(unsigned long) (cigar[i] >> BAM_CIGAR_SHIFT));
		strncat(c,buf,strlen(buf));

		if((cigar[i]&BAM_CIGAR_MASK)>=strlen(CIGAR_TYPES))	// Error
			return 0;

		strncat(c,&(CIGAR_TYPES[cigar[i] & BAM_CIGAR_MASK]),1);
	}
	return strlen(c);
}

static R_INLINE uint8_t *alloc_data(bam1_t *b, int size)
{
	if (b->m_data < size) {
		b->m_data = size;
		kroundup32(b->m_data);
		b->data = (uint8_t*)realloc(b->data, b->m_data);
	}
	return b->data;
}


bam_header_t* clone_bam_header(bam_header_t *h)
{
	bam_header_t *ans=(bam_header_t*)calloc(1, sizeof(bam_header_t));
	ans->n_targets=h->n_targets;
	ans->l_text=h->l_text;
	//ans->n_text=h->n_text;

	ans->text=(char*) calloc(1,(h->l_text)+1);
	strncpy(ans->text,h->text,h->l_text);
	sam_header_parse(ans);
	bam_init_header_hash(ans);
	return ans;
}


SEXP is_nil_externalptr(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[is_nil_externalptr] No external pointer");
		return R_NilValue;
	}
	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=(R_ExternalPtrAddr(ptr)==NULL);
	UNPROTECT(1);
	return ans;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// bam_reader
///////////////////////////////////////////////////////////////////////////////////////////////
static void finalize_bam_reader(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_bam_reader] No external pointer!");
		return;
	}
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(ptr));
	if (reader) {
		samclose(reader);	// checks for 0
		R_SetExternalPtrAddr(ptr,NULL);
	}
}

SEXP bam_reader_open(SEXP filename)
{
	if(TYPEOF(filename)!=STRSXP)
	{
		error("[bam_reader_open] Filename must be a string.\n");
		return R_NilValue;
	}
	const char* _filename=CHAR(STRING_ELT(filename,0));
	samfile_t *reader=samopen(_filename,"rb",0);
	if(!reader)
		error("[bam_reader_open] Opening bam_file \"%s\" failed!",_filename);

	SEXP pReader;
	PROTECT(pReader=R_MakeExternalPtr( (void*)(reader),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(pReader,finalize_bam_reader);
	UNPROTECT(1);
	return pReader;
}

SEXP bam_reader_close(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_close] No external pointer!");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	samclose(reader);
	R_SetExternalPtrAddr(pReader,NULL);
	return R_NilValue;
}


SEXP bam_reader_get_header_text(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_get_header_text] No external pointer!");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	SEXP ans;
	PROTECT(ans=Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(reader->header->text));
	UNPROTECT(1);
	return ans;
}

SEXP bam_reader_get_ref_count(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_get_ref_count] No external pointer!");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=reader->header->n_targets;
	UNPROTECT(1);
	return ans;
}

SEXP bam_reader_get_ref_data(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_get_ref_data] No external pointer!");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_header_t* header=reader->header;

	// create data.frame
	int nProtected=0;
	int nCols=3;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int nRows=header->n_targets;

	// Column 0: ID (RefID)
	SEXP RefID_vector;
	PROTECT(RefID_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: SN (RefName)
	SEXP RefName_vector;
	PROTECT(RefName_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 2: LN (RefLength)
	SEXP RefLength_vector;
	PROTECT(RefLength_vector=allocVector(INTSXP,nRows));
	++nProtected;

	int i;
	for(i=0;i<nRows;++i)
	{
		INTEGER(RefID_vector)[i]=i;
		SET_STRING_ELT(RefName_vector,i,mkChar(header->target_name[i]));
		INTEGER(RefLength_vector)[i]=header->target_len[i];
	}
	SET_VECTOR_ELT(dflist,0,RefID_vector);
	SET_VECTOR_ELT(dflist,1,RefName_vector);
	SET_VECTOR_ELT(dflist,2,RefLength_vector);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;
	SET_STRING_ELT(col_names,0,mkChar("ID"));
	SET_STRING_ELT(col_names,1,mkChar("SN"));
	SET_STRING_ELT(col_names,2,mkChar("LN"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;
    char c[20];
    for(i=0;i<nRows;++i)
    {
    	sprintf(c,"%i",i+1);
    	SET_STRING_ELT(row_names,i,mkChar(c));
    }
    setAttrib(dflist,R_RowNamesSymbol,row_names);
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}

SEXP bam_reader_create_index(SEXP pBamFile,SEXP pIdxFile)
{
	if(TYPEOF(pBamFile)!=STRSXP)
	{
		error("[bam_reader_create_index] BamFile must be a string!\n");
		return R_NilValue;
	}
	if(TYPEOF(pIdxFile)!=STRSXP)
	{
		error("[bam_reader_create_index] IndexFile must be a string!\n");
		return R_NilValue;
	}
	const char *bamFile=CHAR(STRING_ELT(pBamFile,0));
	const char *idxFile=CHAR(STRING_ELT(pIdxFile,0));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=bam_index_build2(bamFile,idxFile);
	UNPROTECT(1);
	return ans;
}

static void finalize_bam_index(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_bam_index] No external pointer!");
		return;
	}
	bam_index_t *index=(bam_index_t *)(R_ExternalPtrAddr(ptr));
	if (index) {
		bam_index_destroy(index);	// checks for zero
		R_SetExternalPtrAddr(ptr,NULL);
	}
}

SEXP bam_reader_load_index(SEXP pIdxFile)
{
	if(TYPEOF(pIdxFile)!=STRSXP)
	{
		error("[bam_reader_load_index] pIdxFile must be a string!\n");
		return R_NilValue;
	}
	const char *idxFile=CHAR(STRING_ELT(pIdxFile,0));
	FILE *f=fopen(idxFile,"rb");
	bam_index_t *index = bam_index_load_core(f);
	fclose(f);

	SEXP idx;
	PROTECT(idx=R_MakeExternalPtr( (void*)(index),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(idx,finalize_bam_index);
	UNPROTECT(1);
	return idx;
}

SEXP bam_reader_unload_index(SEXP pIdx)
{
	if(TYPEOF(pIdx)!=EXTPTRSXP)
	{
		error("[bam_reader_unload_index] No external pointer!\n");
		return R_NilValue;
	}
	bam_index_t *idx=(bam_index_t *)(R_ExternalPtrAddr(pIdx));
	bam_index_destroy(idx);
	R_SetExternalPtrAddr(pIdx,NULL);
	return R_NilValue;
}

SEXP bam_reader_get_next_align(SEXP pReader, SEXP alignedOnly, SEXP primaryOnly)
{
	if(TYPEOF(pReader)!=EXTPTRSXP) {
		error("[bam_reader_get_next_align] No external pointer!\n");
		return R_NilValue;
	}
	if (TYPEOF(alignedOnly) != LGLSXP) {
		error("[bam_reader_get_next_align] alignOnly must be logical!");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam1_t *align=bam_init1();
	_Bool onlyAligned = *LOGICAL( AS_LOGICAL( alignedOnly));
	_Bool onlyPrimary = *LOGICAL( AS_LOGICAL( primaryOnly));
	int naligns = 0;

	while ( naligns < 1) {

		int res=samread(reader,align);
		if(res==-1)
		{
			Rprintf("[getNextAlign] samread found EOF.\n");
			return R_NilValue;
		}
		if(res==-2)
		{
			error("[getNextAlign] samread found truncated BAM-file.\n");
			return R_NilValue;
		}

		if(align==NULL)
			return R_NilValue;
		if (onlyAligned) {
			if ( align->core.flag & BAM_FUNMAP) continue;
		}
		if (onlyPrimary) {
			if ( align->core.flag & BAM_FSECONDARY) continue;
		}
		naligns++;
	}

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)(align),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_align);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_reader_save_aligns(SEXP pReader,SEXP pWriter)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_save_aligns] pReader: No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pWriter)!=EXTPTRSXP)
	{
		error("[bam_reader_save_aligns] pWriter: No external pointer!\n");
		return R_NilValue;
	}

	// reader
	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	if (reader == 0 || !(reader->type & SAM_TYPE_READ))
		error("[bam_reader_save_aligns] Reader not open for reading!");

	// writer
	samfile_t *writer=(samfile_t*)R_ExternalPtrAddr(pWriter);
	if (writer == 0 || (writer->type & SAM_TYPE_READ))
		error("[bam_reader_save_aligns] Writer not open for writing!");


	bam1_t *align=bam_init1();
	unsigned long nAligns=0;
	int res=samread(reader,align);
	while(res>0)
	{
		samwrite(writer,align);
		++nAligns;
		if(nAligns % 1000000 ==0)
			Rprintf("[bam_reader_save_aligns] %u aligns written.\n",nAligns);
		res=samread(reader,align);
	}

	bam_destroy1(align);	// checks for >0!
	if(res==-2)
	{
		error("[bam_reader_save_aligns] samread found truncated BAM-file.\n");
		return R_NilValue;
	}

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=nAligns;
	UNPROTECT(1);
	return ans;
}

SEXP bam_reader_sort_file(SEXP pFilename,SEXP pPrefix,SEXP pMaxMem,SEXP pByName)
{
	if(TYPEOF(pFilename)!=STRSXP)
	{
		error("[bam_writer_sort_file] Filename must be a string\n");
		return R_NilValue;
	}
	if(TYPEOF(pPrefix)!=STRSXP)
	{
		error("[bam_writer_sort_file] Prefix must be a string\n");
		return R_NilValue;
	}
	if(TYPEOF(pMaxMem)!=REALSXP)
	{
		error("[bam_writer_sort_file] MaxMem must be integer value!\n");
		return R_NilValue;
	}
	if(TYPEOF(pByName)!=LGLSXP)
	{
		error("[bam_writer_sort_file] ByName must be bool value!\n");
		return R_NilValue;
	}
	const char *filename=CHAR(STRING_ELT(pFilename,0));
	const char *prefix=CHAR(STRING_ELT(pPrefix,0));
	size_t max_mem=*REAL(pMaxMem);
	_Bool sort_by_name =*(LOGICAL(AS_LOGICAL(pByName)));
	if(sort_by_name)
	{
		bam_sort_core_ext(1, filename, prefix, max_mem, 0);
	}
	else
	{
		bam_sort_core_ext(0, filename, prefix, max_mem, 0);
	}
	return R_NilValue;
}


SEXP bam_reader_get_header(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_get_header] No external pointer!\n");
		return R_NilValue;
	}

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)clone_bam_header(reader->header),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_header);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_reader_tell(SEXP pReader)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_tell] No external pointer!\n");
		return R_NilValue;
	}
	return Rf_ScalarReal(bam_tell(((samfile_t*)(R_ExternalPtrAddr(pReader)))->x.bam));
}

SEXP bam_reader_seek(SEXP pReader, SEXP pPos)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_reader_seek] No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pPos)!=REALSXP)
	{
		error("[bam_reader_seek] Position must be numeric!\n");
		return R_NilValue;
	}
	if(LENGTH(pPos)>1)
	{
		error("[bam_reader_seek] Length of position must be 1!\n");
		return R_NilValue;
	}

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	double *pos=REAL(pPos);
    int64_t res=bam_seek(reader->x.bam,(int64_t)*pos,SEEK_SET);
    if(res<0)
    	Rprintf("[bam_reader_seek] bam_seek fails!\n");
	return R_NilValue;
}


SEXP bam_reader_get_next_chunk(SEXP pReader, SEXP n, SEXP alignedOnly, SEXP primaryOnly)
{
	if(TYPEOF(pReader)!=EXTPTRSXP) {
		error("[bam_reader_get_next_chunk] No external pointer!");
		return R_NilValue;
	}
	if (TYPEOF(n) != INTSXP) {
		error("[bam_reader_get_next_chunk] N must be integer!");
		return R_NilValue;
	}
	if (TYPEOF(alignedOnly) != LGLSXP) {
		error("[bam_reader_get_next_chunk] alignOnly must be logical!");
		return R_NilValue;
	}

	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	int nwant = *INTEGER( AS_INTEGER( n));
	_Bool onlyAligned = *LOGICAL( AS_LOGICAL( alignedOnly));
	_Bool onlyPrimary = *LOGICAL( AS_LOGICAL( primaryOnly));

	align_list *l = init_align_list();
	bam1_t *align = bam_init1();
	int naligns = 0;
	int res;

	while ( naligns < nwant) {
		res = samread( reader, align);
		if( res < 0) break;

		if (onlyAligned) {
			if ( align->core.flag & BAM_FUNMAP) continue;
		}
		if (onlyPrimary) {
			if ( align->core.flag & BAM_FSECONDARY) continue;
		}

		naligns++;
		align_list_push_back( l, align);
	}

	if (onlyAligned && !onlyPrimary && res >= 0) {

		// do a look-ahead to find all the seconary alignments for this read
		bamFile bamfp = reader->x.bam;
		int64_t curpos = bam_tell( bamfp);
		while (1) {
			res = samread( reader, align);
			if( res < 0) break;
			if ( align->core.flag & BAM_FUNMAP) continue;
			if ( align->core.flag & BAM_FSECONDARY) {
				// if this is also a secondary align, keep it too and keep looking
				naligns++;
				align_list_push_back( l, align);
				curpos = bam_tell( bamfp);
			} else {
				// otherwise its a new primary alignment, leave it here for next get
				break;
			}
		}
		// when done with look-ahead, reset to get this same alignment next time
		bam_seek( bamfp, curpos, SEEK_SET);
	}

	// Wind back the result
	l->curr_el=NULL;

    	SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_bam_range);
	UNPROTECT(1);
	return list;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// BamRange
///////////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_range(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_bam_range] No external pointer!\n");
		return;
	}
	align_list *l=(align_list *)(R_ExternalPtrAddr(ptr));
	if (l) {
		destroy_align_list(l);
		l=NULL;
		R_SetExternalPtrAddr(ptr,NULL);
	}
}


SEXP bam_range_init()
{
	align_list *l=init_align_list();
	SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_bam_range);
	UNPROTECT(1);
	return list;
}

static int range_fetch_func(const bam1_t *b, void *data)
{
	align_list *l=(align_list*)data;
	align_list_push_back(l,b);
	return 0;
}

static int range_fetch_complex_func(const bam1_t *b,void *data)
{
	align_list *l=(align_list*)data;
	if(b->core.n_cigar>1)
		align_list_push_back(l,b);
	return 0;
}

SEXP bam_range_fetch(SEXP pReader,SEXP pIndex,SEXP pCoords,SEXP pComplex)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_range_fetch] pReader is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pIndex)!=EXTPTRSXP)
	{
		error("[bam_range_fetch] pIndex is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pCoords)!=REALSXP)
	{
		error("[bam_range_fetch] pCoords is no REAL!\n");
		return R_NilValue;
	}
	if(LENGTH(pCoords)!=3)
	{
		error("[bam_range_fetch] pCoords must contain three values (refid,begin,end)!\n");
		return R_NilValue;
	}
	if(TYPEOF(pComplex)!=LGLSXP)
	{
		error("[bam_range_fetch] pComplex must be logical!\n");
		return R_NilValue;
	}

	samfile_t *reader=(samfile_t*)(R_ExternalPtrAddr(pReader));
	bam_index_t *index=(bam_index_t*)(R_ExternalPtrAddr(pIndex));
	if(reader==NULL)
	{
		error("[bam_range_fetch] Reader must not be NULL pointer!\n");
		return R_NilValue;
	}
	if(index==NULL)
	{
		error("[bam_range_fetch] Index must not be NULL pointer!\n");
		return R_NilValue;
	}

	double *pi=REAL(pCoords);
	int refid=(int) pi[0];
	int begin=(int) pi[1];
	int end=(int) pi[2];

	if(refid<0 || refid >=(reader->header->n_targets))
	{
		error("[bam_range_fetch] refid out of range!\n");
		return R_NilValue;
	}
	if(begin<0 || begin>=end || end>(reader->header->target_len[refid]))
	{
		error("[bam_range_fetch] Begin or end out of range!\n");
		return R_NilValue;
	}

	align_list *l=init_align_list();

	// Retrieve only complex aligns (nCigar>1) when pComplex is set:
	if(LOGICAL(pComplex)[0]==TRUE)
		bam_fetch(reader->x.bam, index, refid, begin, end, (void*)l, range_fetch_complex_func);
	else
		bam_fetch(reader->x.bam, index, refid, begin, end, (void*)l, range_fetch_func);

	// Wind back
    l->curr_el=NULL;

    SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(l),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_bam_range);
	UNPROTECT(1);
	return list;
}

SEXP bam_range_get_next_align(SEXP pRange, SEXP alignedOnly, SEXP primaryOnly)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_next_align] No external pointer!");
	if (TYPEOF(alignedOnly) != LGLSXP) {
		error("[bam_reader_get_next_chunk] alignOnly must be logical!");
		return R_NilValue;
	}
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	const bam1_t *align;
	_Bool onlyAligned = *LOGICAL( AS_LOGICAL( alignedOnly));
	_Bool onlyPrimary = *LOGICAL( AS_LOGICAL( primaryOnly));

	int naligns = 0;

	while ( naligns < 1) {

		align=get_next_align(l);
		if(align==NULL)
			return R_NilValue;
		if (onlyAligned) {
			if ( align->core.flag & BAM_FUNMAP) continue;
		}
		if (onlyPrimary) {
			if ( align->core.flag & BAM_FSECONDARY) continue;
		}
		naligns++;
	}


	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)(align),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_align);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_range_get_prev_align(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_prev_align] No external pointer!\n");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	bam1_t *align=get_prev_align(l);
	if(align==NULL)
		return R_NilValue;

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)(align),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_align);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_range_step_next_align(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_step_next_align] No external pointer!\n");
	pp_curr_align((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}

SEXP bam_range_step_prev_align(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_step_prev_align] No external pointer!\n");
	mm_curr_align((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}


SEXP bam_range_write_current_align(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is No external pointer!\n");
		return R_NilValue;
	}
	write_current_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_insert_past_curr_align(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is No external pointer!\n");
		return R_NilValue;
	}
	insert_past_curr_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_insert_pre_curr_align(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is No external pointer!\n");
		return R_NilValue;
	}
	insert_pre_curr_align((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_mv_curr_align(SEXP pSrc, SEXP pTarget)
{
	// Moves current align in src to end of target list
	// and Moves current align in src to next align
	if(TYPEOF(pSrc)!=EXTPTRSXP)
	{
		error("[bam_range_mv_curr_align] pSrc is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pTarget)!=EXTPTRSXP)
	{
		error("[bam_range_mv_curr_align] pTarget is No external pointer!\n");
		return R_NilValue;
	}

	align_list_mv_curr_elem((align_list*)(R_ExternalPtrAddr(pSrc)),(align_list*)(R_ExternalPtrAddr(pTarget)));
	return R_NilValue;
}


SEXP bam_range_get_align_df(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_align_df] No external pointer!");

	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	// Read Adress of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);
	const bam1_t *align;

	// create data.frame
	int nProtected=0;
	int nCols=10;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int nRows=(l->size);
	int i,j;

	// Column 0: refid
	SEXP ref_vector;
	PROTECT(ref_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: position
	SEXP pos_vector;
	PROTECT(pos_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 2: nCigar
	SEXP nCigar_vector;
	PROTECT(nCigar_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 3: cigar
	SEXP cig_vector;
	PROTECT(cig_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 4: flag
	SEXP flag_vector;
	PROTECT(flag_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 5: seq
	SEXP seq_vector;
	PROTECT(seq_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 6: qual
	SEXP qual_vector;
	PROTECT(qual_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 7: name
	SEXP name_vector;
	PROTECT(name_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 8: strand_reverse
	SEXP strev_vector;
	PROTECT(strev_vector=allocVector(LGLSXP,nRows));
	++nProtected;

	// Column 9: insert size
	SEXP isize_vector;
	PROTECT(isize_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// seq+cigar
	unsigned char *raw_seq;
	int32_t seq_len;
	int buf_size=2048;
	char *buf=(char*) calloc(buf_size,sizeof(char));
	uint8_t *quals;

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		INTEGER(ref_vector)[i]=(align->core.tid);
		INTEGER(pos_vector)[i]=(align->core.pos) + BASE_1;
		INTEGER(isize_vector)[i]=(align->core.isize);
		INTEGER(nCigar_vector)[i]=(align->core.n_cigar);

		/////////////////////////////////////////
		// Cigar String
		if(cigar2str(buf,align)==0)
		{
			error("[bam_align_get_align_df] Cigar error!\n");
			return R_NilValue;
		}
		SET_STRING_ELT(cig_vector,i,mkChar(buf));
		//Rprintf("%i\t%s\n",i,buf);
		clear_buf(buf,buf_size);
		/////////////////////////////////////////

		INTEGER(flag_vector)[i]=(align->core.flag);
		/////////////////////////////////////////
		// seq
		seq_len=align->core.l_qseq;
		if(seq_len>buf_size)
		{
			buf_size=2*(seq_len+1);
			free(buf);
			buf= (char*) calloc(buf_size,sizeof(char));
		}
		raw_seq=bam1_seq(align);
		for(j=0;j<seq_len;++j)
			buf[j]=bam_nt16_rev_table[bam1_seqi(raw_seq,j)];
		buf[j]=0;

		SET_STRING_ELT(seq_vector,i,mkChar(buf));
		////////////////////////////////////////
		// quals
			quals=bam1_qual(align);
			for(j=0;j<seq_len;++j)
				buf[j]=(char) (quals[j]+33);
			buf[j]=0;
			SET_STRING_ELT(qual_vector,i,mkChar(buf));

		/////////////////////////////////////////
		// read name
		SET_STRING_ELT(name_vector,i,mkChar(bam1_qname(align)));

		////////////////////////////////////////
		// strand_info
		LOGICAL(strev_vector)[i]=bam1_strand(align);
	}

	// Reset curr_el pointer
	l->curr_el=e;

	SET_VECTOR_ELT(dflist,0,ref_vector);
	SET_VECTOR_ELT(dflist,1,pos_vector);
	SET_VECTOR_ELT(dflist,2,nCigar_vector);
	SET_VECTOR_ELT(dflist,3,cig_vector);
	SET_VECTOR_ELT(dflist,4,flag_vector);
	SET_VECTOR_ELT(dflist,5,seq_vector);
	SET_VECTOR_ELT(dflist,6,qual_vector);
	SET_VECTOR_ELT(dflist,7,name_vector);
	SET_VECTOR_ELT(dflist,8,strev_vector);
	SET_VECTOR_ELT(dflist,9,isize_vector);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("refid"));
	SET_STRING_ELT(col_names,1,mkChar("position"));
	SET_STRING_ELT(col_names,2,mkChar("nCigar"));
	SET_STRING_ELT(col_names,3,mkChar("cigar"));
	SET_STRING_ELT(col_names,4,mkChar("flag"));
	SET_STRING_ELT(col_names,5,mkChar("seq"));
	SET_STRING_ELT(col_names,6,mkChar("qual"));
	SET_STRING_ELT(col_names,7,mkChar("name"));
	SET_STRING_ELT(col_names,8,mkChar("revstrand"));
	SET_STRING_ELT(col_names,9,mkChar("insertsize"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;

    for(i=0;i<nRows;++i)
    {
    	sprintf(buf,"%i",i+1);
    	SET_STRING_ELT(row_names,i,mkChar(buf));
    }
    free(buf);
    setAttrib(dflist,R_RowNamesSymbol,row_names);
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}

SEXP bam_range_write(SEXP pWriter,SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_write] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pWriter)!=EXTPTRSXP)
	{
		error("[bam_range_write] pWriter is No external pointer!\n");
		return R_NilValue;
	}

	unsigned long bytes_written=0;
	unsigned long range_size, i;

	samfile_t *writer=(samfile_t*) R_ExternalPtrAddr(pWriter);
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	range_size=l->size;
	// For restoration
	align_element *e=l->curr_el;
	wind_back(l);

	// Retrieve const pointer (i.e. no copy)
	for(i=0;i<range_size;++i)
		bytes_written+=samwrite(writer,get_const_next_align(l));

	// restore curr_el
	l->curr_el=e;

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=bytes_written;
	UNPROTECT(1);
	return ans;
}

SEXP bam_range_wind_back(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_wind_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	wind_back(l);
	return R_NilValue;
}

SEXP bam_range_get_size(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_get_size] pRange is No external pointer!\n");
		return R_NilValue;
	}
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=(l->size);
	UNPROTECT(1);
	return ans;
}

SEXP bam_range_push_back(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_back] pAlign is No external pointer!\n");
		return R_NilValue;
	}
	align_list_push_back((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_push_front(SEXP pRange,SEXP pAlign)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_push_front] pRange is No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_range_push_front] pAlign is No external pointer!\n");
		return R_NilValue;
	}
	align_list_push_front((align_list*)(R_ExternalPtrAddr(pRange)),(bam1_t*) (R_ExternalPtrAddr(pAlign)));
	return R_NilValue;
}

SEXP bam_range_pop_back(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_pop_back] pRange is No external pointer!\n");
		return R_NilValue;
	}
	align_list_pop_back((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}
SEXP bam_range_pop_front(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
	{
		error("[bam_range_pop_front] pRange is No external pointer!\n");
		return R_NilValue;
	}
	align_list_pop_front((align_list*)(R_ExternalPtrAddr(pRange)));
	return R_NilValue;
}

SEXP bam_range_get_readID(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_readID] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	int nRows=(l->size);
	const bam1_t *align;
	int i;

	SEXP readIDans;
	PROTECT(readIDans=allocVector(STRSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		SET_STRING_ELT(readIDans,i,mkChar(bam1_qname(align)));
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	return readIDans;
}


SEXP bam_range_get_refid(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_refid] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	int nRows=(l->size);
	const bam1_t *align;
	int i;

	SEXP refIDans;
	PROTECT(refIDans=allocVector(INTSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		INTEGER(refIDans)[i]=align->core.tid;
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	return refIDans;
}


SEXP bam_range_get_position(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_position] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	int nRows=(l->size);
	const bam1_t *align;
	int i;

	SEXP positionans;
	PROTECT(positionans=allocVector(INTSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		INTEGER(positionans)[i]=(align->core.pos) + BASE_1;
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	return positionans;
}


SEXP bam_range_get_insert_size(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_insert_size] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	int nRows=(l->size);
	const bam1_t *align;
	int i;

	SEXP sizeans;
	PROTECT(sizeans=allocVector(INTSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		INTEGER(sizeans)[i]=(align->core.isize);
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	return sizeans;
}


SEXP bam_range_get_tag(SEXP pRange, SEXP tag)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_tag] No external pointer!");
	if(TYPEOF(tag)!=STRSXP)
		error("[bam_range_get_tag] tag must be a string!!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	const char *tg = CHAR( STRING_ELT( tag, 0));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	int buf_len = 2048;
	char *buf = (char*) calloc( buf_len, sizeof(char));
	int nRows=(l->size);
	int i;

	SEXP tagans;
	PROTECT(tagans=allocVector(STRSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		strcpy( buf, bam_aux_format1( bam_aux_get_core( align, tg)));
		if (strlen(buf) >= buf_len) error( "[bam_align_get_tag] Buffer overflow!");
		SET_STRING_ELT(tagans,i,mkChar(buf));
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	free(buf);
	return tagans;
}

SEXP bam_range_get_all_tags(SEXP pRange, SEXP sep)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_all_tags] No external pointer!");
	if(TYPEOF(sep)!=STRSXP)
		error("[bam_range_get_all_tags] sep must be a string!!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	const char *separ = CHAR( STRING_ELT( sep, 0));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	int buf_len = 4096;
	char *buf = (char*) calloc( buf_len, sizeof(char));
	int nRows=(l->size);
	int i;

	SEXP tagans;
	PROTECT(tagans=allocVector(STRSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		strcpy( buf, bam_aux_formatAll( align, separ));
		if (strlen(buf) >= buf_len) error( "[bam_align_get_all_tags] Buffer overflow!");
		SET_STRING_ELT(tagans,i,mkChar(buf));
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	free(buf);
	return tagans;
}

SEXP bam_range_set_tag(SEXP pRange, SEXP tag, SEXP value)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_set_tag] No external pointer!");
	if(TYPEOF(tag)!=STRSXP)
		error("[bam_range_set_tag] tag must be a string!!");
	if(TYPEOF(value)!=STRSXP)
		error("[bam_range_set_tag] value must be a string!!");

	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	const char *tg = CHAR( STRING_ELT( tag, 0));
	const char *val;
	const bam1_t *align;

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	int nRows=(l->size);
	int nVal = LENGTH( value);
	int i, lval;
	char buf[1024];

	if ( nRows != nVal) 
		error( "[bam_range_set_tag] length of values not equal to size of range!");

	//FILE *fp = fopen( "setTag.txt", "w");
	//fprintf( fp, "tag = %s \n", tg);

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		//align=get_next_align(l);
		val = CHAR( STRING_ELT( value, i));
		lval = strlen(val);
		memcpy( buf, val, lval);
		buf[lval] = 0;

		// first see if this tag is already here
		uint8_t *oldtag = bam_aux_get( align, tg);

		if ( oldtag) {
			//fprintf( fp, "I: %d  Found: %s   Set: %s \n", i, (char*)oldtag, buf);
			bam_aux_del( (bam1_t*) align, oldtag);
//		} else {
//			fprintf( fp, "I: %d  Set: %s \n", i, buf);
		}
//		fflush(fp);

		// OK, now add this new tag, with the trailing null char
		if (lval) bam_aux_append( (bam1_t*) align, tg, 'Z', lval+1, (uint8_t*) buf);

		//write_current_align( l, align);
	}
	//fclose( fp);
	// Reset curr_el pointer
	l->curr_el=e;

	return R_NilValue;
}


SEXP bam_range_get_read_sequence(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_read_seq] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	unsigned char *raw_seq;
	int buf_len = 4096;
	char *seq= (char*) calloc(buf_len,sizeof(char));
	char *seq2= (char*) calloc(buf_len,sizeof(char));
	int nRows=(l->size);
	int32_t i,j,seq_len;

	SEXP seqans;
	PROTECT(seqans=allocVector(STRSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		raw_seq=bam1_seq(align);
		seq_len=align->core.l_qseq;
		if (seq_len >= buf_len) error( "[bam_align_get_read_sequence] Buffer overflow!");
		for(j=0;j<seq_len;++j)
			seq[j]=bam_nt16_rev_table[bam1_seqi(raw_seq,j)];
		seq[j]=0;
		if ( bam1_strand(align)) {
			strcpy( seq2, seq);
			reversecomplement( seq, seq2);
		}
		SET_STRING_ELT(seqans,i,mkChar(seq));
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	free(seq2);
	free(seq);
	return seqans;
}

SEXP bam_range_get_align_sequence(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_align_seq] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	unsigned char *raw_seq;
	int buf_len = 4096;
	char *seq= (char*) calloc(buf_len,sizeof(char));
	int nRows=(l->size);
	int32_t i,j,seq_len;

	SEXP seqans;
	PROTECT(seqans=allocVector(STRSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		raw_seq=bam1_seq(align);
		seq_len=align->core.l_qseq;
		if (seq_len >= buf_len) error( "[bam_align_get_align_sequence] Buffer overflow!");
		for(j=0;j<seq_len;++j)
			seq[j]=bam_nt16_rev_table[bam1_seqi(raw_seq,j)];
		seq[j]=0;
		SET_STRING_ELT(seqans,i,mkChar(seq));
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	free(seq);
	return seqans;
}

SEXP bam_range_get_align_qualities(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_align_qualities] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	unsigned char *raw_qual;
	int buf_len = 4096;
	char *qual= (char*) calloc(buf_len,sizeof(char));
	int nRows=(l->size);
	int32_t i,j,qual_len;

	SEXP qualans;
	PROTECT(qualans=allocVector(STRSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		raw_qual=bam1_qual(align);
		qual_len=align->core.l_qseq;
		if (qual_len >= buf_len) error( "[bam_align_get_align_qualities] Buffer overflow!");
		for(j=0;j<qual_len;++j)
			qual[j]=raw_qual[j] + 33;
		qual[j]=0;
		SET_STRING_ELT(qualans,i,mkChar(qual));
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	free(qual);
	return qualans;
}

SEXP bam_range_get_read_qualities(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_read_qualities] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	unsigned char *raw_qual;
	int buf_len = 4096;
	char *qual= (char*) calloc(buf_len,sizeof(char));
	char *qual2= (char*) calloc(buf_len,sizeof(char));
	int nRows=(l->size);
	int32_t i,j,qual_len;

	SEXP qualans;
	PROTECT(qualans=allocVector(STRSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		raw_qual=bam1_qual(align);
		qual_len=align->core.l_qseq;
		if (qual_len >= buf_len) error( "[bam_align_get_read_qualities] Buffer overflow!");
		for(j=0;j<qual_len;++j)
			qual[j]=raw_qual[j] + 33;
		qual[j]=0;
		if ( bam1_strand(align)) {
			strcpy( qual2, qual);
			reversequality( qual, qual2);
		}

		SET_STRING_ELT(qualans,i,mkChar(qual));
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	free(qual2);
	free(qual);
	return qualans;
}

SEXP bam_range_get_read_phredScores(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_read_qualities] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	unsigned char *raw_qual;
	int buf_len = 4096;
	char *qual= (char*) calloc(buf_len,sizeof(char));
	char *qual2= (char*) calloc(buf_len,sizeof(char));
	int nRows=(l->size);
	int32_t i,j,qual_len,this_qual_len;

	// use the first to get the read length
	align=get_const_next_align(l);
	raw_qual=bam1_qual(align);
	qual_len=align->core.l_qseq;
	wind_back(l);

	SEXP qualans;
	PROTECT(qualans=allocMatrix(INTSXP,nRows,qual_len));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		raw_qual=bam1_qual(align);
		this_qual_len=align->core.l_qseq;
		if (this_qual_len >= buf_len) error( "[bam_align_get_read_phredScores] Buffer overflow!");
		if (this_qual_len > qual_len) this_qual_len <- qual_len;
		for(j=0;j<this_qual_len;++j)
			qual[j]=raw_qual[j];
		qual[j]=0;
		if ( bam1_strand(align)) {
			strcpy( qual2, qual);
			reversequality( qual, qual2);
		}

		for(j=0;j<this_qual_len;++j)
			INTEGER(qualans)[ j*nRows + i] = (int) qual[j];
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	free(qual2);
	free(qual);
	return qualans;
}

SEXP bam_range_get_map_quality(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_map_quality] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	int nRows=(l->size);
	int32_t i;

	SEXP mapans;
	PROTECT(mapans=Rf_allocVector(INTSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		INTEGER(mapans)[i]=align->core.qual;
	}
	UNPROTECT(1);

	// Reset curr_el pointer
	l->curr_el=e;

	return mapans;
}


SEXP bam_range_get_align_length(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_align_length] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	int nRows=(l->size);
	int32_t i;

	SEXP mapans;
	PROTECT(mapans=Rf_allocVector(INTSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		INTEGER(mapans)[i]=align->core.l_qseq;
	}
	UNPROTECT(1);

	// Reset curr_el pointer
	l->curr_el=e;

	return mapans;
}


SEXP bam_range_get_strand_reverse(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP) {
		error("[bam_range_strand_reverse] No external pointer!");
		return R_NilValue;
	}
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	int nRows=(l->size);
	int32_t i;

	SEXP boolans;
	PROTECT(boolans=allocVector(LGLSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		LOGICAL(boolans)[i]= bam1_strand( align);
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	return boolans;
}


SEXP bam_range_is_first_in_pair(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_is_first_in_pair] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	int nRows=(l->size);
	int32_t i;

	SEXP boolans;
	PROTECT(boolans=allocVector(LGLSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		LOGICAL(boolans)[i]= ((align->core.flag & BAM_FREAD1) != 0);
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	return boolans;
}

SEXP bam_range_is_second_in_pair(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_is_second_in_pair] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	int nRows=(l->size);
	int32_t i;

	SEXP boolans;
	PROTECT(boolans=allocVector(LGLSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		LOGICAL(boolans)[i]=  ((align->core.flag & BAM_FREAD2) != 0);
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	return boolans;
}

SEXP bam_range_get_read_locus(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_read_locus] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one data frame
	int nRows=(l->size);
	int nCols=4;
	SEXP dframe;
	PROTECT(dframe=allocVector(VECSXP,nCols));
	int nProtect=1;

	int32_t i;
	const bam1_t *align;
	int values[4];

	SEXP lane;
	PROTECT(lane=allocVector(INTSXP, nRows));
	SEXP tile;
	PROTECT(tile=allocVector(INTSXP, nRows));
	SEXP x;
	PROTECT(x=allocVector(INTSXP, nRows));
	SEXP y;
	PROTECT(y=allocVector(INTSXP, nRows));
	nProtect += 4;

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		read_locus( bam1_qname(align), (int*) &values);
		INTEGER(lane)[i] = values[0];
		INTEGER(tile)[i] = values[1];
		INTEGER(x)[i] = values[2];
		INTEGER(y)[i] = values[3];
	}
	// Reset curr_el pointer
	l->curr_el=e;

	SET_VECTOR_ELT(dframe, 0, lane);
	SET_VECTOR_ELT(dframe, 1, tile);
	SET_VECTOR_ELT(dframe, 2, x);
	SET_VECTOR_ELT(dframe, 3, y);

	SEXP names;
	PROTECT( names=allocVector(STRSXP,4));
	nProtect++;
	SET_STRING_ELT( names,0,mkChar("lane"));
	SET_STRING_ELT( names,1,mkChar("tile"));
	SET_STRING_ELT( names,2,mkChar("x"));
	SET_STRING_ELT( names,3,mkChar("y"));
	setAttrib( dframe, R_NamesSymbol, names);

	SEXP rownames;
	PROTECT( rownames=allocVector(STRSXP,nRows));
	nProtect++;
	char *buf = calloc( 32, sizeof(char));
	for( i=0; i<nRows; i++) {
		sprintf( buf, "%i", (i+1));
		SET_STRING_ELT(rownames, i, mkChar(buf));
	}
	free(buf);
	setAttrib( dframe, R_RowNamesSymbol, rownames);
	setAttrib( dframe, R_ClassSymbol, mkString("data.frame"));
	UNPROTECT(nProtect);
	return dframe;
}

SEXP bam_range_is_unmapped(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_is_unmapped] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	int nRows=(l->size);
	int32_t i;

	SEXP boolans;
	PROTECT(boolans=allocVector(LGLSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		LOGICAL(boolans)[i]= ((align->core.flag & BAM_FUNMAP) != 0);
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	return boolans;
}

SEXP bam_range_is_secondary_align(SEXP pRange)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_is_unmapped] No external pointer!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create one vector
	const bam1_t *align;
	int nRows=(l->size);
	int32_t i;

	SEXP boolans;
	PROTECT(boolans=allocVector(LGLSXP,nRows));

	for(i=0;i<nRows;++i)
	{
		align=get_const_next_align(l);
		LOGICAL(boolans)[i]= ((align->core.flag & BAM_FSECONDARY) != 0);
	}
	// Reset curr_el pointer
	l->curr_el=e;

	UNPROTECT(1);
	return boolans;
}

SEXP bam_range_get_cigar_df(SEXP pRange, SEXP readUnits)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_cigar_df] No external pointer!");
	if(TYPEOF(readUnits)!=LGLSXP)
		error("[bam_range_get_cigar_df] 'readUnits' logical!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	_Bool useReadUnits = *LOGICAL( AS_LOGICAL( readUnits));

	const bam1_t *align;
	int nRows=(l->size);

	// create data.frame
	int nProtected=0;
	int nCols=4;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	int i,j,k;

	// allocate the columns, guess a max returned length
	int nOut = nRows * MAX_CIGAR_PER_ALIGN;
	SEXP Pos_vector;
	PROTECT(Pos_vector=allocVector(INTSXP,nOut));
	++nProtected;
	SEXP Length_vector;
	PROTECT(Length_vector=allocVector(INTSXP,nOut));
	++nProtected;
	SEXP ID_vector;
	PROTECT(ID_vector=allocVector(INTSXP,nOut));
	++nProtected;
	SEXP Type_vector;
	PROTECT(Type_vector=allocVector(STRSXP,nOut));
	++nProtected;

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	int nout=0;
	for(k=0;k<nRows;++k)
	{
		align=get_const_next_align(l);
		int nRows=align->core.n_cigar;
		uint32_t *cigar=bam1_cigar(align);
		int curPos = 0;
		for(i=0;i<nRows;++i)
		{
			j = i;
			if (useReadUnits && bam1_strand(align))  j = (nRows-1-i);
			if((cigar[j]&BAM_CIGAR_MASK)>=strlen(CIGAR_TYPES))
			{
				error("[bam_range_getCigar_df] Cigar_type not in defined range!");
				return R_NilValue;
			}
			int thisLen=cigar[j] >> BAM_CIGAR_SHIFT;
			INTEGER(Pos_vector)[nout]=curPos + BASE_1;
			INTEGER(Length_vector)[nout]=thisLen;
			curPos += thisLen;
			SET_STRING_ELT(Type_vector,nout,mkCharLen(CIGAR_TYPES+(cigar[j]&BAM_CIGAR_MASK),1));
			INTEGER(ID_vector)[nout]=k+1;
			nout++;
		}
	}
	SETLENGTH( ID_vector, nout);
	SETLENGTH( Pos_vector, nout);
	SETLENGTH( Length_vector, nout);
	SETLENGTH( Type_vector, nout);

	// Reset curr_el pointer
	l->curr_el=e;

	SET_VECTOR_ELT(dflist,0,ID_vector);
	SET_VECTOR_ELT(dflist,1,Type_vector);
	SET_VECTOR_ELT(dflist,2,Pos_vector);
	SET_VECTOR_ELT(dflist,3,Length_vector);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("Align"));
	SET_STRING_ELT(col_names,1,mkChar("Type"));
	SET_STRING_ELT(col_names,2,mkChar("Position"));
	SET_STRING_ELT(col_names,3,mkChar("Length"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	SEXP row_names;
	PROTECT(row_names=allocVector(STRSXP,nout));
	++nProtected;

	char c[20];
	for(i=0;i<nout;++i)
	{
		sprintf(c,"%i",i+1);
		SET_STRING_ELT(row_names,i,mkChar(c));
	}
	setAttrib(dflist,R_RowNamesSymbol,row_names);

	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}

SEXP bam_range_modify(SEXP pRange, SEXP RefID, SEXP Pos, 
				SEXP ReadID, SEXP Seq, SEXP Qual)
{
	if(TYPEOF(pRange)!=EXTPTRSXP) {
		error("[bam_range_modify] No external pointer!\n");
		return R_NilValue;
	}
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	align_list *newl = init_align_list();
	const bam1_t *oldAlign;
	bam1_t *newAlign;


	// now see if we were passed any variable pieces that will need to be repacked
	int haveNew = 0;
	if ( TYPEOF(ReadID) == STRSXP) {
		haveNew=1;
	}
	if ( TYPEOF(Seq) == STRSXP) {
		haveNew=1;
	}
	if ( TYPEOF(Qual) == STRSXP) {
		haveNew=1;
	}

	int nRows=(l->size);
	int32_t i;

	for(i=0;i<nRows;++i)
	{
		oldAlign=get_const_next_align(l);
		newAlign=bam_dup1( oldAlign);

		// fixed fields that are not NULL
		if ( TYPEOF(RefID) == INTSXP) {
			int newRefID = INTEGER( RefID)[i];
			newAlign->core.tid=newRefID;
		}
		if ( TYPEOF(Pos) == INTSXP) {
			int newPos = INTEGER( Pos)[i];
			newAlign->core.pos=newPos - BASE_1;
		}


		if ( haveNew) {
			int lqname=oldAlign->core.l_qname;
			char qname[1024];
			memcpy( qname, bam1_qname(oldAlign), lqname);
			int lcigar=oldAlign->core.n_cigar;
			int lcigar_bytes = lcigar*4;
			char cigar[1024];
			memcpy( cigar, bam1_cigar(oldAlign), lcigar_bytes);
			int lqseq=oldAlign->core.l_qseq;
			int lqseq_bytes = (lqseq+1)/2;
			char qseq[1024];
			memcpy( qseq, bam1_seq(oldAlign), lqseq_bytes);
			int lqual=oldAlign->core.l_qseq;
			char qual[1024];
			memcpy( qual, bam1_qual(oldAlign), lqseq);
			int laux=oldAlign->l_aux;
			char aux[4096];
			memcpy( aux, bam1_aux(oldAlign), laux);
	
			// now see what we were passed
			if ( TYPEOF(ReadID) == STRSXP) {
				strcpy( qname, CHAR( STRING_ELT(ReadID,i)));
				lqname=strlen(qname) + 1;
			}
			if ( TYPEOF(Seq) == STRSXP) {
				char qseqtmp[1024];
				strcpy( qseqtmp, CHAR( STRING_ELT(Seq,i)));
				lqseq=strlen(qseqtmp);
				memset( qseq, 0, lqseq);
				for( int j=0; j<lqseq; j++) 
					qseq[j/2] |= bam_nt16_table[ (int)qseqtmp[j]] << 4*(1-j%2);
				lqseq_bytes=(lqseq+1)/2;
			}
			if ( TYPEOF(Qual) == STRSXP) {
				strcpy( qual, CHAR( STRING_ELT(Qual,i)));
				lqual=strlen(qual);
				for (int j=0;j<lqual;j++) {
					qual[j] = qual[j] - 33;
				}
			}
			if (lqual != lqseq) {
				error( "[bam_align_modify]:  Length of SEQ and QUAL strings must be equal!");
				return R_NilValue;
			}
	
			// build the new variable data packet
			// force the CIGAR string to be compatible... call it all M for now...
			if ( lqseq != oldAlign->core.l_qseq) {
				newAlign->core.n_cigar = lcigar = 1;
				uint32_t tmp_cigar = lqseq << BAM_CIGAR_SHIFT | BAM_CMATCH;
				lcigar_bytes = lcigar * 4;
				memcpy( cigar, &tmp_cigar, lcigar_bytes);
			}
			int newDataLen = lqname + lcigar_bytes + lqseq_bytes + lqseq + laux;
			newAlign->data = (uint8_t*) alloc_data( newAlign, newDataLen);
			char *p = newAlign->data;
			memcpy( p, qname, lqname);
			p += lqname;
			memcpy( p, cigar, lcigar_bytes);
			p += lcigar_bytes;
			memcpy( p, qseq, lqseq_bytes);
			p += lqseq_bytes;
			memcpy( p, qual, lqseq);
			p += lqseq;
			memcpy( p, aux, laux);
			// update those count fields
			newAlign->core.l_qname=lqname;
			newAlign->core.l_qseq=lqseq;
			newAlign->data_len=newDataLen;
		}
		align_list_push_back( newl, newAlign);
	}

	// Wind back
	l->curr_el=NULL;

    	SEXP list;
	PROTECT(list=R_MakeExternalPtr( (void*)(newl),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(list,finalize_bam_range);
	UNPROTECT(1);
	return list;
}
	

SEXP bam_range_get_mismatch_df(SEXP pRange, SEXP readUnits)
{
	if(TYPEOF(pRange)!=EXTPTRSXP)
		error("[bam_range_get_mismatch_df] No external pointer!");
	if(TYPEOF(readUnits)!=LGLSXP)
		error("[bam_range_get_mismatch_df] 'readUnits' logical!");
	align_list *l=(align_list*)(R_ExternalPtrAddr(pRange));
	_Bool useReadUnits = *LOGICAL( AS_LOGICAL( readUnits));

	// Read Address of current align->for reconstitution at the end
	align_element *e=l->curr_el;
	wind_back(l);

	// create data.frame
	int nProtected=0;
	int nCols=5;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	int nRows=(l->size);
	// estimate a size based on CIGAR & MISMATCH counts
	int nOut = nRows * MAX_CIGAR_PER_ALIGN * MAX_MISMATCH_PER_ALIGN / 8;

	// Column 0,1,2: Type, Position & Length
	SEXP Type_vector;
	PROTECT(Type_vector=allocVector(STRSXP,nOut));
	++nProtected;
	SEXP Pos_vector;
	PROTECT(Pos_vector=allocVector(INTSXP,nOut));
	++nProtected;
	SEXP ReadBase_vector;
	PROTECT(ReadBase_vector=allocVector(STRSXP,nOut));
	++nProtected;
	SEXP RefBase_vector;
	PROTECT(RefBase_vector=allocVector(STRSXP,nOut));
	++nProtected;
	SEXP ID_vector;
	PROTECT(ID_vector=allocVector(INTSXP,nOut));
	++nProtected;

	int i,j,k;
	const bam1_t *align;
	char seq[1024];
	char mmbuf[4096];
	char refbase[512], readbase[512];
	int thisCIGlen, curCIGpos, thisMMlen, curMMpos, outpos;
	char thisCIGchar, thisMMchar, thisType;
	int needCIG, needMM;
	int nCigar, curPos, mmlength, seq_len, icig, loopcnt;
	uint32_t *cigar;
	unsigned char *raw_seq;
	_Bool doRC;
	char *curMM;

	// loop over all alignments
	int nout = 0;
	for(k=0;k<nRows;++k) {

	align=get_const_next_align(l);

	// we need to track both the CIGAR
	nCigar=align->core.n_cigar;

	cigar=bam1_cigar(align);
	curPos = 0;
	// and the MD tag at the same time
	strcpy( mmbuf, bam_aux_format1( bam_aux_get( align, BAMTAG_MISMATCH)));
	mmlength = strlen( mmbuf);
	// and the aligned sequence
	seq_len = align->core.l_qseq;
	raw_seq=bam1_seq(align);
	for(i=0;i<seq_len;++i)
		seq[i]=bam_nt16_rev_table[bam1_seqi(raw_seq,i)];
	seq[i]=0;
	doRC = (useReadUnits && bam1_strand(align));

	// ready to run along all at once...
	curMM = mmbuf;
	curCIGpos = curMMpos = 0;
	icig = 0;
	// get the first entry from both CIGAR and MM
	needCIG = needMM = 1;
	loopcnt = 0;
	while ( icig < nCigar || curMM <= (mmbuf+mmlength)) {
		//printf( "\nDebug WHILE:  %d %d %d %d %d %c %c %d %d ", nout, needCIG, needMM, icig, 
		//		(curMM-mmbuf), thisCIGchar, thisMMchar, curCIGpos, curMMpos);
		if ( loopcnt++ > 100) break;
		if (needCIG) {
			if((cigar[icig]&BAM_CIGAR_MASK)>=strlen(CIGAR_TYPES)) {
				error("[bam_align_getCigar_df] Cigar_type not in defined range!");
				return R_NilValue;
			}
			thisCIGlen = cigar[icig] >> BAM_CIGAR_SHIFT;
			strncpy( &thisCIGchar, (CIGAR_TYPES+(cigar[icig]&BAM_CIGAR_MASK)), 1);
			curCIGpos += thisCIGlen;
			icig++;
			//printf( "\ngot CIG: (nout, thisCIGlen, thisCIGchar, curCIGpos)  %d  %d  %c  %d ", nout, 
			//		thisCIGlen, thisCIGchar, curCIGpos);
			needCIG=0;
		}
		if (needMM) {
			thisMMlen = (int) strtol( curMM, &curMM, 10);
			curMMpos += thisMMlen;
			thisMMchar = *curMM++;
			//printf( "\ngot MM: (nout, thisMMlen, thisMMchar, curMMpos)  %d  %d  %c  %d ", nout, thisMMlen, 
			//		thisMMchar, curMMpos);
			needMM=0;
		}
		// a deletion should be mentioned in both worlds
		if (thisMMchar == '^' && thisCIGchar == 'D') {
			thisType = 'D';
			for ( int ii=0; ii<thisCIGlen; ii++) refbase[ii] = *curMM++;
			refbase[thisCIGlen] = 0;
			readbase[0] = 0;
			outpos = curMMpos;
			// the deletion was from the reference, so don't step ahead on this position tracker
			curCIGpos -= thisCIGlen;
			//printf( "\n Got 'Delete' (curCIGpos, curMMpos, refbase, readbase, outpos)  %d  %d  %s %s %d ", curCIGpos, 
			//		curMMpos, refbase, readbase, outpos);
			needMM = 1;
			SET_STRING_ELT( Type_vector, nout, mkCharLen( &thisType, 1));
			SET_STRING_ELT( RefBase_vector, nout, mkChar( refbase));
			SET_STRING_ELT( ReadBase_vector, nout, mkChar( readbase));
			INTEGER( Pos_vector)[nout] = outpos + BASE_1;
			if (doRC) INTEGER( Pos_vector)[nout] = (seq_len - outpos);
			INTEGER( ID_vector)[nout] = k+1;
			nout++;
			continue;
		}
		// the default thing in the mismatch string is mismatches
		if (curMMpos < curCIGpos) {
			thisType = 'X';
			refbase[0] = thisMMchar;
			refbase[1] = 0;
			readbase[0] = seq[curMMpos];
			readbase[1] = 0;
			if (doRC) reversecomplement( readbase, readbase);
			outpos = curMMpos;
			//printf( "\n Got 'Mismatch' (curCIGpos, curMMpos, ref, read, outpos)  %d  %d  %s  %s  %d ", curCIGpos, 
			//		curMMpos, refbase, readbase, outpos);
			needMM = 1;
			// this mismatch covered one base
			curMMpos++;
			SET_STRING_ELT( Type_vector, nout, mkCharLen( &thisType, 1));
			SET_STRING_ELT( RefBase_vector, nout, mkChar( refbase));
			SET_STRING_ELT( ReadBase_vector, nout, mkChar( readbase));
			INTEGER( Pos_vector)[nout] = outpos + BASE_1;
			if (doRC) INTEGER( Pos_vector)[nout] = (seq_len - outpos);
			INTEGER( ID_vector)[nout] = k+1;
			nout++;
			continue;
		}
		// the insertion flag is ONLY in the CIGAR
		if ( thisCIGchar == 'I') {
			thisType = 'I';
			outpos = curCIGpos - thisCIGlen; // it starts 'before' its extent
			for ( int ii=0; ii<thisCIGlen; ii++) readbase[ii] = seq[outpos+ii];
			readbase[thisCIGlen] = 0;
			refbase[0] = 0;
			//printf( "\n Got 'Insert' (curCIGpos, curMMpos, refbase, readbase, outpos)  %d  %d  %s %s %d ", curCIGpos, 
			//		curMMpos, refbase, readbase, outpos);
			needCIG = 1;
			curMMpos += thisCIGlen;
			SET_STRING_ELT( Type_vector, nout, mkCharLen( &thisType, 1));
			SET_STRING_ELT( RefBase_vector, nout, mkChar( refbase));
			SET_STRING_ELT( ReadBase_vector, nout, mkChar( readbase));
			INTEGER( Pos_vector)[nout] = outpos + BASE_1;
			if (doRC) INTEGER( Pos_vector)[nout] = (seq_len - outpos);
			INTEGER( ID_vector)[nout] = k+1;
			nout++;
			continue;
		}
		// no explicit change, see who we need to step ahead on...
		if ( curMMpos >= curCIGpos) {
			needCIG = 1;
			continue;
		}
	} // end of each single alignment record
	} // end of all aligns in this chunk

	// Wind back
	l->curr_el=NULL;


	SETLENGTH( ID_vector, nout);
	SETLENGTH( Type_vector, nout);
	SETLENGTH( Pos_vector, nout);
	SETLENGTH( RefBase_vector, nout);
	SETLENGTH( ReadBase_vector, nout);

	SET_VECTOR_ELT(dflist,0,ID_vector);
	SET_VECTOR_ELT(dflist,1,Type_vector);
	SET_VECTOR_ELT(dflist,2,Pos_vector);
	SET_VECTOR_ELT(dflist,3,RefBase_vector);
	SET_VECTOR_ELT(dflist,4,ReadBase_vector);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("Align"));
	SET_STRING_ELT(col_names,1,mkChar("Type"));
	SET_STRING_ELT(col_names,2,mkChar("Position"));
	SET_STRING_ELT(col_names,3,mkChar("RefBase"));
	SET_STRING_ELT(col_names,4,mkChar("ReadBase"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	SEXP row_names;
	PROTECT(row_names=allocVector(STRSXP,nout));
	++nProtected;

	char c[20];
	for(i=0;i<nout;++i)
	{
		sprintf(c,"%i",i+1);
		SET_STRING_ELT(row_names,i,mkChar(c));
	}
	setAttrib(dflist,R_RowNamesSymbol,row_names);

	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}


///////////////////////////////////////////////////////////////////////////////////////////////
// bam_header
///////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_header(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_bam_header] No external pointer!");
		return;
	}
	bam_header_t* header=(bam_header_t*)(R_ExternalPtrAddr(ptr));
	if(header) {
		bam_header_destroy(header);
		R_SetExternalPtrAddr(ptr,NULL);
	}
}

SEXP init_bam_header(SEXP pHeaderText)
{
	if(TYPEOF(pHeaderText)==NILSXP)
		return R_NilValue;

	if(TYPEOF(pHeaderText)!=STRSXP)
	{
		error("[init_bam_header] Header Text must be a string.\n");
		return R_NilValue;
	}
	bam_header_t *h=(bam_header_t*)calloc(1, sizeof(bam_header_t));

	// Copy header text
	const char* header_text=CHAR(STRING_ELT(pHeaderText,0));
	h->l_text=strlen(header_text);
	char *txt=(char*) calloc(1,(h->l_text)+1);
	strncpy(txt,CHAR(STRING_ELT(pHeaderText,0)),h->l_text);
	h->text=txt;

	sam_header_parse(h);
	bam_init_header_hash(h);

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)h,R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_header);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_header_get_header_text(SEXP pHeader)
{
	if(TYPEOF(pHeader)!=EXTPTRSXP)
	{
		error("[bam_header_get_header_text] No external pointer!");
		return R_NilValue;
	}

	bam_header_t* header=(bam_header_t*)(R_ExternalPtrAddr(pHeader));
	SEXP ans;
	PROTECT(ans=Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(header->text));
	UNPROTECT(1);
	return ans;
}
///////////////////////////////////////////////////////////////////////////////////////////////
// bam_writer
///////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_writer(SEXP ptr)
{
	if(TYPEOF(ptr)!=EXTPTRSXP)
	{
		error("[finalize_bam_writer] No external pointer!");
		return;
	}
	samfile_t *writer=(samfile_t*)(R_ExternalPtrAddr(ptr));
	if(writer)
	{
		samclose(writer);
		R_SetExternalPtrAddr(ptr,NULL);
		Rprintf("[bamWriter] finalized.\n");
	}
}

SEXP bam_writer_open(SEXP pHeader,SEXP pFilename)
{
	if(TYPEOF(pHeader)!=EXTPTRSXP)
	{
		error("[bam_writer_open] pHeader No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pFilename)!=STRSXP)
	{
		error("[bam_writer_open] pFilename no string!\n");
		return R_NilValue;
	}

	bam_header_t *header=(bam_header_t*) (R_ExternalPtrAddr(pHeader));
	samfile_t    *writer=samopen(CHAR(STRING_ELT(pFilename,0)),"wb",header);
	if(!writer)
		error("[bam_writer_open] samopen returned NULL pointer!\n");

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr( (void*) (writer),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_writer);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_reader_open_writer(SEXP pReader,SEXP pFilename)
{
	if(TYPEOF(pReader)!=EXTPTRSXP)
	{
		error("[bam_writer_open] pReader No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pFilename)!=STRSXP)
	{
		error("[bam_writer_open] pFilename no string!\n");
		return R_NilValue;
	}
	samfile_t *reader=(samfile_t*) (R_ExternalPtrAddr(pReader));
	samfile_t *writer=samopen(CHAR(STRING_ELT(pFilename,0)),"wb",reader->header);

	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr( (void*) (writer),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_writer);
	UNPROTECT(1);
	return ptr;
}

SEXP bam_writer_save_align(SEXP pWriter, SEXP pAlign)
{
	if(TYPEOF(pWriter)!=EXTPTRSXP)
	{
		error("[bam_writer_save_align] No external pointer!\n");
		return R_NilValue;
	}
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_writer_save_align] No external pointer!\n");
		return R_NilValue;
	}
	samfile_t *writer=(samfile_t*)R_ExternalPtrAddr(pWriter);
	bam1_t *align=(bam1_t*)R_ExternalPtrAddr(pAlign);

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=samwrite(writer,align);
	UNPROTECT(1);
	return ans;
}

SEXP bam_writer_close(SEXP pWriter)
{
	if(TYPEOF(pWriter)!=EXTPTRSXP)
	{
		error("[bam_writer_close] No exteranl pointer!\n");
		return R_NilValue;
	}
	samfile_t *writer= (samfile_t*) (R_ExternalPtrAddr(pWriter));
	samclose(writer);
	R_SetExternalPtrAddr(pWriter,NULL);
	Rprintf("[bamWriter] closed.\n");
	return R_NilValue;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// bam_align
///////////////////////////////////////////////////////////////////////////////////////////////

static void finalize_bam_align(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[finalize_bam_align] No external pointer!");
	bam1_t *align= (bam1_t*)(R_ExternalPtrAddr(pAlign));
	if (align) {
		bam_destroy1(align);	// checks for >0!
		R_SetExternalPtrAddr(pAlign,NULL);
	}
}


SEXP bam_align_modify(SEXP pAlign, SEXP RefID, SEXP Pos, 
				SEXP ReadID, SEXP Seq, SEXP Qual)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP) {
		error("[bam_align_modify] No external pointer!\n");
		return R_NilValue;
	}

	bam1_t *oldAlign=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	bam1_t *newAlign=bam_dup1( oldAlign);

	// fixed fields that are not NULL
	if ( TYPEOF(RefID) == INTSXP || TYPEOF(RefID) == REALSXP) {
		int newRefID = *INTEGER( AS_INTEGER( RefID));
		newAlign->core.tid=newRefID;
	}
	if ( TYPEOF(Pos) == INTSXP || TYPEOF(RefID) == REALSXP) {
		int newPos = *INTEGER( AS_INTEGER( Pos));
		newAlign->core.pos=newPos - BASE_1;
	}

	// now see if we were passed any variable pieces that will need to be repacked
	int haveNew = 0;
	if ( TYPEOF(ReadID) == STRSXP) {
		haveNew=1;
	}
	if ( TYPEOF(Seq) == STRSXP) {
		haveNew=1;
	}
	if ( TYPEOF(Qual) == STRSXP) {
		haveNew=1;
	}

	if ( haveNew) {
		int lqname=oldAlign->core.l_qname;
		char qname[1024];
		memcpy( qname, bam1_qname(oldAlign), lqname);
		int lcigar=oldAlign->core.n_cigar;
		int lcigar_bytes = lcigar*4;
		char cigar[1024];
		memcpy( cigar, bam1_cigar(oldAlign), lcigar_bytes);
		int lqseq=oldAlign->core.l_qseq;
		int lqseq_bytes = (lqseq+1)/2;
		char qseq[1024];
		memcpy( qseq, bam1_seq(oldAlign), lqseq_bytes);
		int lqual=oldAlign->core.l_qseq;
		char qual[1024];
		memcpy( qual, bam1_qual(oldAlign), lqseq);
		int laux=oldAlign->l_aux;
		char aux[4096];
		memcpy( aux, bam1_aux(oldAlign), laux);

		// now see what we were passed
		if ( TYPEOF(ReadID) == STRSXP) {
			strcpy( qname, CHAR( STRING_ELT(ReadID,0)));
			lqname=strlen(qname) + 1;
		}
		if ( TYPEOF(Seq) == STRSXP) {
			char qseqtmp[1024];
			strcpy( qseqtmp, CHAR( STRING_ELT(Seq,0)));
			lqseq=strlen(qseqtmp);
			memset( qseq, 0, lqseq);
			for( int j=0; j<lqseq; j++) 
				qseq[j/2] |= bam_nt16_table[ (int)qseqtmp[j]] << 4*(1-j%2);
			lqseq_bytes=(lqseq+1)/2;
		}
		if ( TYPEOF(Qual) == STRSXP) {
			strcpy( qual, CHAR( STRING_ELT(Qual,0)));
			lqual=strlen(qual);
			for (int j=0;j<lqual;j++) {
				qual[j] = qual[j] - 33;
			}
		}
		if (lqual != lqseq) {
			error( "[bam_align_modify]:  Length of SEQ and QUAL strings must be equal!");
			return R_NilValue;
		}

		// build the new variable data packet
		// force the CIGAR string to be compatible... call it all M for now...
		if ( lqseq != oldAlign->core.l_qseq) {
			newAlign->core.n_cigar = lcigar = 1;
			uint32_t tmp_cigar = lqseq << BAM_CIGAR_SHIFT | BAM_CMATCH;
			lcigar_bytes = lcigar * 4;
			memcpy( cigar, &tmp_cigar, lcigar_bytes);
		}
		int newDataLen = lqname + lcigar_bytes + lqseq_bytes + lqseq + laux;
		newAlign->data = (uint8_t*) alloc_data( newAlign, newDataLen);
		char *p = newAlign->data;
		memcpy( p, qname, lqname);
		p += lqname;
		memcpy( p, cigar, lcigar_bytes);
		p += lcigar_bytes;
		memcpy( p, qseq, lqseq_bytes);
		p += lqseq_bytes;
		memcpy( p, qual, lqseq);
		p += lqseq;
		memcpy( p, aux, laux);
		// update those count fields
		newAlign->core.l_qname=lqname;
		newAlign->core.l_qseq=lqseq;
		newAlign->data_len=newDataLen;
	}

	// send it back
	SEXP ptr;
	PROTECT(ptr=R_MakeExternalPtr((void*)(newAlign),R_NilValue,R_NilValue));
	R_RegisterCFinalizer(ptr,finalize_bam_align);
	UNPROTECT(1);
	return ptr;
}


SEXP bam_align_get_readID(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_readID] No external pointer!");
	bam1_t *align= (bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(bam1_qname(align)));
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_read_locus(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_read_locus] No external pointer!");
	bam1_t *align= (bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,4));

	int values[4];
	read_locus( bam1_qname(align), (int*) &values);
	for( int i=0; i<4; i++) INTEGER(ans)[i] = values[i];

	SEXP names;
	PROTECT( names=allocVector(STRSXP,4));
	SET_STRING_ELT( names,0,mkChar("lane"));
	SET_STRING_ELT( names,1,mkChar("tile"));
	SET_STRING_ELT( names,2,mkChar("x"));
	SET_STRING_ELT( names,3,mkChar("y"));
	setAttrib( ans, R_NamesSymbol, names);

	UNPROTECT(2);
	return ans;
}

SEXP bam_align_get_refid(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_getRefID] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.tid;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_position(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_position] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=(align->core.pos) + BASE_1;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_nCigar(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_nCigar] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.n_cigar;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_cigar_df(SEXP pAlign, SEXP readUnits)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_cigar_df] No external pointer!");
	if(TYPEOF(readUnits)!=LGLSXP)
		error("[bam_align_get_cigar_df] 'readUnits' logical!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	_Bool useReadUnits = *LOGICAL( AS_LOGICAL( readUnits));

	// create data.frame
	int nProtected=0;
	int nCols=3;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int nRows=align->core.n_cigar;
	int i,j;

	// Column 1,2: Position & Length
	SEXP Pos_vector;
	PROTECT(Pos_vector=allocVector(INTSXP,nRows));
	++nProtected;
	SEXP Length_vector;
	PROTECT(Length_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 0: Type
	SEXP Type_vector;
	PROTECT(Type_vector=allocVector(STRSXP,nRows));
	++nProtected;

	uint32_t *cigar=bam1_cigar(align);
	int curPos = 0;
	for(i=0;i<nRows;++i)
	{
		j = i;
		if (useReadUnits && bam1_strand(align))  j = (nRows-1-i);
		if((cigar[j]&BAM_CIGAR_MASK)>=strlen(CIGAR_TYPES))
		{
			error("[bam_align_getCigar_df] Cigar_type not in defined range!");
			return R_NilValue;
		}
		int thisLen=cigar[j] >> BAM_CIGAR_SHIFT;
		INTEGER(Pos_vector)[i]=curPos + BASE_1;
		INTEGER(Length_vector)[i]=thisLen;
		curPos += thisLen;
		SET_STRING_ELT(Type_vector,i,mkCharLen(CIGAR_TYPES+(cigar[j]&BAM_CIGAR_MASK),1));
	}

	SET_VECTOR_ELT(dflist,0,Type_vector);
	SET_VECTOR_ELT(dflist,1,Pos_vector);
	SET_VECTOR_ELT(dflist,2,Length_vector);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("Type"));
	SET_STRING_ELT(col_names,1,mkChar("Position"));
	SET_STRING_ELT(col_names,2,mkChar("Length"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	SEXP row_names;
	PROTECT(row_names=allocVector(STRSXP,nRows));
	++nProtected;

	char c[20];
	for(i=0;i<nRows;++i)
	{
		sprintf(c,"%i",i+1);
		SET_STRING_ELT(row_names,i,mkChar(c));
	}
	setAttrib(dflist,R_RowNamesSymbol,row_names);

	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}

SEXP bam_align_get_mismatch_df(SEXP pAlign, SEXP readUnits)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_mismatch_df] No external pointer!");
	if(TYPEOF(readUnits)!=LGLSXP)
		error("[bam_align_get_mismatch_df] 'readUnits' logical!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	_Bool useReadUnits = *LOGICAL( AS_LOGICAL( readUnits));

	// create data.frame
	int nProtected=0;
	int nCols=4;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	int i,j;

	int nCigar=align->core.n_cigar;

	// Column 0,1,2: Type, Position & Length
	int nOut = nCigar * MAX_MISMATCH_PER_ALIGN;
	SEXP Type_vector;
	PROTECT(Type_vector=allocVector(STRSXP,nOut));
	++nProtected;
	SEXP Pos_vector;
	PROTECT(Pos_vector=allocVector(INTSXP,nOut));
	++nProtected;
	SEXP ReadBase_vector;
	PROTECT(ReadBase_vector=allocVector(STRSXP,nOut));
	++nProtected;
	SEXP RefBase_vector;
	PROTECT(RefBase_vector=allocVector(STRSXP,nOut));
	++nProtected;

	// we need to track both the CIGAR
	uint32_t *cigar=bam1_cigar(align);
	int curPos = 0;
	// and the MD tag at the same time
	char mmbuf[1024];
	strcpy( mmbuf, bam_aux_format1( bam_aux_get( align, BAMTAG_MISMATCH)));
	int mmlength = strlen( mmbuf);
	// and the aligned sequence
	int seq_len = align->core.l_qseq;
	char *seq= (char*) calloc(seq_len+1,sizeof(char));
	unsigned char *raw_seq=bam1_seq(align);
	for(i=0;i<seq_len;++i)
		seq[i]=bam_nt16_rev_table[bam1_seqi(raw_seq,i)];
	seq[i]=0;
	//char *rcseq= (char*) calloc(seq_len+1,sizeof(char));
	//reversecomplement( rcseq, seq);
	_Bool doRC = (useReadUnits && bam1_strand(align));
	//if (doRC) strcpy( seq, rcseq);
	char refbase[128], readbase[128];

	// ready to run along all at once...
	int thisCIGlen, curCIGpos, thisMMlen, curMMpos, outpos;
	char thisCIGchar, thisMMchar, thisType;
	char *curMM = mmbuf;
	curCIGpos = curMMpos = 0;
	int icig = 0;
	int nout = 0;
	// get the first entry from both CIGAR and MM
	int needCIG, needMM;
	needCIG = needMM = 1;
	int loopcnt = 0;
	int DEBUG = 0;
	while ( icig < nCigar || curMM <= (mmbuf+mmlength)) {
		if (DEBUG) printf( "\nDebug WHILE:  %d %d %d %d %d %c %c %d %d ", nout, needCIG, needMM, icig, 
				(curMM-mmbuf), thisCIGchar, thisMMchar, curCIGpos, curMMpos);
		if ( loopcnt++ > 100) break;
		if (needCIG) {
			if((cigar[icig]&BAM_CIGAR_MASK)>=strlen(CIGAR_TYPES)) {
				error("[bam_align_getCigar_df] Cigar_type not in defined range!");
				return R_NilValue;
			}
			thisCIGlen = cigar[icig] >> BAM_CIGAR_SHIFT;
			strncpy( &thisCIGchar, (CIGAR_TYPES+(cigar[icig]&BAM_CIGAR_MASK)), 1);
			curCIGpos += thisCIGlen;
			icig++;
			if (DEBUG) printf( "\ngot CIG: (nout, thisCIGlen, thisCIGchar, curCIGpos)  %d  %d  %c  %d ", nout, 
					thisCIGlen, thisCIGchar, curCIGpos);
			needCIG=0;
		}
		if (needMM) {
			thisMMlen = (int) strtol( curMM, &curMM, 10);
			curMMpos += thisMMlen;
			thisMMchar = *curMM++;
			if (DEBUG) printf( "\ngot MM: (nout, thisMMlen, thisMMchar, curMMpos)  %d  %d  %c  %d ", nout, thisMMlen, 
					thisMMchar, curMMpos);
			needMM=0;
		}
		// a deletion should be mentioned in both worlds
		if (thisMMchar == '^' && thisCIGchar == 'D') {
			thisType = 'D';
			for ( int ii=0; ii<thisCIGlen; ii++) refbase[ii] = *curMM++;
			refbase[thisCIGlen] = 0;
			readbase[0] = 0;
			outpos = curMMpos;
			// the deletion was from the reference, so don't step ahead on this position tracker
			curCIGpos -= thisCIGlen;
			if (DEBUG) printf( "\n Got 'Delete' (curCIGpos, curMMpos, refbase, readbase, outpos)  %d  %d  %s %s %d ", curCIGpos, 
					curMMpos, refbase, readbase, outpos);
			needMM = needCIG = 1;
			SET_STRING_ELT( Type_vector, nout, mkCharLen( &thisType, 1));
			SET_STRING_ELT( RefBase_vector, nout, mkChar( refbase));
			SET_STRING_ELT( ReadBase_vector, nout, mkChar( readbase));
			INTEGER( Pos_vector)[nout] = outpos + BASE_1;
			if (doRC) INTEGER( Pos_vector)[nout] = (seq_len - outpos);
			nout++;
			continue;
		}
		// the default thing in the mismatch string is mismatches
		if (curMMpos < curCIGpos) {
			thisType = 'X';
			refbase[0] = thisMMchar;
			refbase[1] = 0;
			readbase[0] = seq[curMMpos];
			readbase[1] = 0;
			if (doRC) reversecomplement( readbase, readbase);
			outpos = curMMpos;
			if (DEBUG) printf( "\n Got 'Mismatch' (curCIGpos, curMMpos, ref, read, outpos)  %d  %d  %s  %s  %d ", curCIGpos, 
					curMMpos, refbase, readbase, outpos);
			needMM = 1;
			// this mismatch covered one base
			curMMpos++;
			SET_STRING_ELT( Type_vector, nout, mkCharLen( &thisType, 1));
			SET_STRING_ELT( RefBase_vector, nout, mkChar( refbase));
			SET_STRING_ELT( ReadBase_vector, nout, mkChar( readbase));
			INTEGER( Pos_vector)[nout] = outpos + BASE_1;
			if (doRC) INTEGER( Pos_vector)[nout] = (seq_len - outpos);
			nout++;
			continue;
		}
		// the insertion flag is ONLY in the CIGAR
		if ( thisCIGchar == 'I') {
			thisType = 'I';
			outpos = curCIGpos - thisCIGlen; // it starts 'before' its extent
			for ( int ii=0; ii<thisCIGlen; ii++) readbase[ii] = seq[outpos+ii];
			readbase[thisCIGlen] = 0;
			refbase[0] = 0;
			if (DEBUG) printf( "\n Got 'Insert' (curCIGpos, curMMpos, refbase, readbase, outpos)  %d  %d  %s %s %d ", curCIGpos, 
					curMMpos, refbase, readbase, outpos);
			needCIG = 1;
			curMMpos += thisCIGlen;
			SET_STRING_ELT( Type_vector, nout, mkCharLen( &thisType, 1));
			SET_STRING_ELT( RefBase_vector, nout, mkChar( refbase));
			SET_STRING_ELT( ReadBase_vector, nout, mkChar( readbase));
			INTEGER( Pos_vector)[nout] = outpos + BASE_1;
			if (doRC) INTEGER( Pos_vector)[nout] = (seq_len - outpos);
			nout++;
			continue;
		}
		// no explicit change, see who we need to step ahead on...
		if ( curMMpos >= curCIGpos) {
			needCIG = 1;
			continue;
		}
	}
	free(seq);

	SETLENGTH( Type_vector, nout);
	SETLENGTH( Pos_vector, nout);
	SETLENGTH( RefBase_vector, nout);
	SETLENGTH( ReadBase_vector, nout);

	SET_VECTOR_ELT(dflist,0,Type_vector);
	SET_VECTOR_ELT(dflist,1,Pos_vector);
	SET_VECTOR_ELT(dflist,2,RefBase_vector);
	SET_VECTOR_ELT(dflist,3,ReadBase_vector);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("Type"));
	SET_STRING_ELT(col_names,1,mkChar("Position"));
	SET_STRING_ELT(col_names,2,mkChar("RefBase"));
	SET_STRING_ELT(col_names,3,mkChar("ReadBase"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	SEXP row_names;
	PROTECT(row_names=allocVector(STRSXP,nout));
	++nProtected;

	char c[20];
	for(i=0;i<nout;++i)
	{
		sprintf(c,"%i",i+1);
		SET_STRING_ELT(row_names,i,mkChar(c));
	}
	setAttrib(dflist,R_RowNamesSymbol,row_names);

	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	return dflist;
}

SEXP bam_align_get_mate_refid(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_mate_refid] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.mtid;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_mate_position(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_mate_position] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=(align->core.mpos) + BASE_1;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_insert_size(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_insert_size] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.isize;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_map_quality(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_map_quality] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.qual;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_align_length(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_align_length] No external pointer!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.l_qseq;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_align_sequence(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_get_segment_sequence] No external pointer!");
		return R_NilValue;
	}

	// Extract char* sequence with samtools
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	int32_t seq_len=align->core.l_qseq;
	char *seq= (char*) calloc(seq_len+1,sizeof(char));
	unsigned char *raw_seq=bam1_seq(align);
	int32_t i;
	for(i=0;i<seq_len;++i)
		seq[i]=bam_nt16_rev_table[bam1_seqi(raw_seq,i)];
	seq[i]=0;

	SEXP ans;
	PROTECT(ans=allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(seq));
	free(seq);
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_read_sequence(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_get_segment_sequence] No external pointer!");
		return R_NilValue;
	}

	// Extract char* sequence with samtools
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	int32_t seq_len=align->core.l_qseq;
	char *seq= (char*) calloc(seq_len+1,sizeof(char));
	char *seq2= (char*) calloc(seq_len+1,sizeof(char));
	unsigned char *raw_seq=bam1_seq(align);
	int32_t i;
	for(i=0;i<seq_len;++i)
		seq[i]=bam_nt16_rev_table[bam1_seqi(raw_seq,i)];
	seq[i]=0;

	if ( bam1_strand(align)) {
		strcpy( seq2, seq);
		reversecomplement( seq, seq2);
	}

	SEXP ans;
	PROTECT(ans=allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(seq));
	UNPROTECT(1);
	free(seq2);
	free(seq);
	return ans;
}

SEXP bam_align_get_align_qualities(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_get_qualities] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	int32_t qual_len=align->core.l_qseq;
	char *qual= (char*) calloc(qual_len+1,sizeof(char));
	char *raw_qual=bam1_qual(align);
	for( int i=0; i<qual_len; i++) {
		qual[i] = (char) (raw_qual[i] + 33);
	}

	SEXP ans;
	PROTECT(ans=allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar( qual));
	UNPROTECT(1);
	free( qual);
	return ans;
}

SEXP bam_align_get_read_qualities(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_get_qualities] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	int32_t qual_len=align->core.l_qseq;
	char *qual= (char*) calloc(qual_len+1,sizeof(char));
	char *qual2= (char*) calloc(qual_len+1,sizeof(char));
	char *raw_qual=bam1_qual(align);
	for( int i=0; i<qual_len; i++) {
		qual[i] = (char) (raw_qual[i] + 33);
	}

	if ( bam1_strand(align)) {
		strcpy( qual2, qual);
		reversequality( qual, qual2);
	}

	SEXP ans;
	PROTECT(ans=allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar(qual));
	UNPROTECT(1);
	free(qual2);
	free(qual);
	return ans;
}

SEXP bam_align_get_read_phredScores(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_get_qualities] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	int32_t qual_len=align->core.l_qseq;
	char *qual= (char*) calloc(qual_len+1,sizeof(char));
	char *qual2= (char*) calloc(qual_len+1,sizeof(char));
	char *raw_qual=bam1_qual(align);
	for( int i=0; i<qual_len; i++) {
		qual[i] = (char) (raw_qual[i]);
	}

	if ( bam1_strand(align)) {
		strcpy( qual2, qual);
		reversequality( qual, qual2);
	}

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,qual_len));
	for( int i=0; i<qual_len; i++) {
		INTEGER(ans)[i] = (int) qual[i];
	}
	UNPROTECT(1);
	free(qual2);
	free(qual);
	return ans;
}

SEXP bam_align_get_tag(SEXP pAlign, SEXP tag)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_tag] No external pointer!");
	if(TYPEOF(tag)!=STRSXP)
		error("[bam_align_get_tag] tag must be a string!!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	const char *tg = CHAR( STRING_ELT( tag, 0));

	char *buf = (char*) calloc( 1024, sizeof(char));
	strcpy( buf, bam_aux_format1( bam_aux_get( align, tg)));

	SEXP ans;
	PROTECT(ans=Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar( buf));
	UNPROTECT(1);
	free( buf);
	return ans;
}


SEXP bam_align_get_all_tags(SEXP pAlign, SEXP sep)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_tag] No external pointer!");
	if(TYPEOF(sep)!=STRSXP)
		error("[bam_align_get_tag] sep must be a string!!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	const char *separator = CHAR( STRING_ELT( sep, 0));

	char *buf = (char*) calloc( 4096, sizeof(char));
	strcpy( buf, bam_aux_formatAll( align, separator));

	uint8_t *s = bam1_aux(align);
	uint8_t *end = align->data + align->data_len;

	SEXP ans;
	PROTECT(ans=Rf_allocVector(STRSXP,1));
	SET_STRING_ELT(ans,0,mkChar( buf));
	UNPROTECT(1);
	free( buf);
	return ans;
}


SEXP bam_align_set_tag(SEXP pAlign, SEXP tag, SEXP value)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
		error("[bam_align_get_tag] No external pointer!");
	if(TYPEOF(tag)!=STRSXP)
		error("[bam_align_get_tag] tag must be a string!!");
	if(TYPEOF(value)!=STRSXP)
		error("[bam_align_get_tag] value must be a string!!");
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	const char *tg = CHAR( STRING_ELT( tag, 0));
	const char *val = CHAR( STRING_ELT( value, 0));
	char buf[1024];
	int lval = strlen(val);
	memcpy( buf, val, lval);
	buf[lval] = 0;

	// first see if this tag is already here
	uint8_t *oldtag = bam_aux_get( align, tg);
	if ( oldtag) bam_aux_del( align, oldtag);

	// OK, now add this new tag, with the trailing null char
	if (lval) bam_aux_append( align, tg, 'Z', lval+1, (uint8_t*) buf);
	return R_NilValue;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Alignment flags

SEXP bam_align_is_paired(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_paired] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align= (bam1_t*) (R_ExternalPtrAddr(pAlign));
	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FPAIRED;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_mapped_in_proper_pair(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_mapped_in_proper_pair] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*) (R_ExternalPtrAddr(pAlign));
	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FPROPER_PAIR;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_is_unmapped(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_unmapped] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FUNMAP;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_mate_is_unmapped(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_mate_is_unmapped] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FUNMAP;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_strand_reverse(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_strand_reverse] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*) (R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=bam1_strand(align);
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_mate_strand_reverse(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_strand_reverse] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=bam1_mstrand(align);
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_is_first_in_pair(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_first_in_pair] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]= ((align->core.flag & BAM_FREAD1) != 0);
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_is_second_in_pair(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_second_in_pair] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]= ((align->core.flag & BAM_FREAD2) != 0);
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_is_secondary_align(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_secondary_align] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FSECONDARY;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_fail_qc(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_fail_qc] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FQCFAIL;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_is_pcr_or_optical_dup(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_is_pcr_or_optical_dup] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(LGLSXP,1));
	LOGICAL(ans)[0]=align->core.flag & BAM_FDUP;
	UNPROTECT(1);
	return ans;
}

SEXP bam_align_get_flag(SEXP pAlign)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_get_flag] No external pointer!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));

	SEXP ans;
	PROTECT(ans=allocVector(INTSXP,1));
	INTEGER(ans)[0]=align->core.flag;
	UNPROTECT(1);
	return ans;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Writing accessors
SEXP bam_align_set_is_paired(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_paired] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_paired] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FPAIRED);
	return R_NilValue;
}

SEXP bam_align_set_mapped_in_proper_pair(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_mapped_in_proper_pair] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_mapped_in_proper_pair] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FPROPER_PAIR);
	return R_NilValue;
}

SEXP bam_align_set_is_unmapped(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_unmapped] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_unmapped] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FUNMAP);
	return R_NilValue;
}

SEXP bam_align_set_mate_is_unmapped(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_mate_is_unmapped] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_mate_is_unmapped] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FMUNMAP);
	return R_NilValue;
}

SEXP bam_align_set_strand_reverse(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_strand_reverse] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_strand_reverse] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),16);
	return R_NilValue;
}

SEXP bam_align_set_mate_strand_reverse(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_mate_strand_reverse] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_mate_strand_reverse] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),32);
	return R_NilValue;
}

SEXP bam_align_set_is_first_in_pair(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_first_in_pair] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_first_in_pair] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FREAD1);
	return R_NilValue;
}

SEXP bam_align_set_is_second_in_pair(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_second_in_pair] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_second_in_pair] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FREAD2);
	return R_NilValue;
}

SEXP bam_align_set_is_secondary_align(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_secondary_align] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_secondary_align] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FSECONDARY);
	return R_NilValue;
}

SEXP bam_align_set_fail_qc(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_fail_qc] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_fail_qc] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FQCFAIL);
	return R_NilValue;
}

SEXP bam_align_set_is_pcr_or_optical_dup(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_is_pcr_or_optical_dup] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=LGLSXP)
	{
		error("[bam_align_set_is_pcr_or_optical_dup] No bool value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	set_flag(align,*(LOGICAL(AS_LOGICAL(val))),BAM_FDUP);
	return R_NilValue;
}

SEXP bam_align_set_flag(SEXP pAlign, SEXP val)
{
	if(TYPEOF(pAlign)!=EXTPTRSXP)
	{
		error("[bam_align_set_flag] No external pointer!");
		return R_NilValue;
	}
	if(TYPEOF(val)!=INTSXP)
	{
		error("[bam_align_set_flag] No integer value!");
		return R_NilValue;
	}
	bam1_t *align=(bam1_t*)(R_ExternalPtrAddr(pAlign));
	align->core.flag=*INTEGER(val);
	return R_NilValue;
}

////////////////////////
// 	extra functions
///////////////////////

void reversequality( char *dest, const char* src) {

	int len = strlen( src);
	int i, j;
	j = len;
	for ( i=0; i<len; i++) {
		dest[ --j] = src[i];
	}
	dest[len] = '\0';
}


void reversecomplement( char *dest, const char* src) {

	int len = strlen( src);
	char *tmp = (char*) calloc( len+2, 1);
	strcpy( tmp, src);
	int i, j;
	char c;
	j = len;
	for ( i=0; i<len; i++) {
		c = dest[ --j] = tmp[i];
		if ( c == 'C') dest[j] = 'G';
		else if ( c == 'G') dest[j] = 'C';
		else if ( c == 'T') dest[j] = 'A';
		else if ( c == 'A') dest[j] = 'T';
		else ;
	}
	free( tmp);
	dest[len] = '\0';
}

void read_locus( const char *readID, int *values)
{
	//  disect the lane:tile:x:y from the end of the readID
	// and write 4 integers into the storage location given
	int sep = (int) ':';
	int val[4] = { 0, 0, 0, 0};
	int i;
	char *p;

	// assume error condition - preload output
	for(i=0; i<4; i++) values[i] = val[i];
	if ( readID == 0) return;
	int str_len = strlen(readID);
	if ( str_len < 7) return;

	char *s = (char*) calloc( str_len + 2, sizeof(char));
	strcpy( s, readID);

	// there 'may' be blank space after the location fields
	p = strchr( s, ' ');
	if (p) *p = '\0';

	// now grab y, x, tile, lane from tail of the ID
	for( i=3; i >= 0; i--) {
		p = strrchr( s, sep);
		if (p) {
			*p = '\0';
			val[i] = atoi( p+1);
		}
	}

	// if all are still zero, try another separator
	if ( (val[0]+val[1]+val[2]+val[3]) == 0) {
		sep = (int) '_';
		for( i=3; i >= 0; i--) {
			p = strrchr( s, sep);
			if (p) {
				*p = '\0';
				val[i] = atoi( p+1);
			}
		}
		// there is a tiny chance that the lane number is the very first character of the readID
		if (val[0] == 0) val[0] = atoi(s);
	}

	// pass them back
	for(i=0; i<4; i++) values[i] = val[i];
}
#endif
