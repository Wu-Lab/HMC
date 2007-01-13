
#include "HaploFile.h"

#include "MemLeak.h"


////////////////////////////////
//
// class HaploFile

HaploFile::HaploFile()
{
	m_filename[0] = 0;
	m_haplo_data = NULL;
	m_has_id = true;
}

HaploFile::HaploFile(const char *filename)
{
	setFileName(filename);
	m_haplo_data = NULL;
	m_has_id = true;
}

void HaploFile::readHaploData(HaploData &hd)
{
	FILE *fp;
	char line[STR_LEN_FILE_LINE], buf[STR_LEN_GENOTYPE_ID];
	char *s, *delim = " \t\r\n";
	Haplotype h1, h2;
	int i, j;
	m_haplo_data = &hd;
	fp = fopen(m_filename, "r");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename);
		exit(1);
	}
	fscanf(fp, "%d\n", &i);
	fscanf(fp, "%d\n", &j);
	if (i <= 0 || j <= 0) {
		Logger::error("Invalid file type!");
		exit(1);
	}
	hd.setGenotypeNum(i);
	hd.setGenotypeLen(j);
	// set loci positions
	fgets(line, STR_LEN_FILE_LINE, fp);
	s = line + strspn(line, delim);
	if (s[0] == 'P') {
		s += strcspn(s, delim);
		s += strspn(s, delim);
		for (i=0; i<hd.genotype_len(); i++) {
			hd.setAllelePosition(i, atoi(s));
			s += strcspn(s, delim);
			s += strspn(s, delim);
		}
		fgets(line, STR_LEN_FILE_LINE, fp);
		s = line + strspn(line, delim);
	}
	// set loci type
	for (i=0; i<hd.genotype_len(); i++) {
		hd.setAlleleType(i, s[0]);
		s++;
		s += strspn(s, delim);
	}
	for (i=0; i<hd.genotype_num(); i++) {
		if (m_has_id) {
			fgets(line, STR_LEN_FILE_LINE, fp);
			if (sscanf(line, "%s", buf) > 0) {
				hd[i].setID(buf);
			}
		}
		fgets(line, STR_LEN_FILE_LINE, fp);
		h1.read(hd.allele_type().c_str(), line, hd.genotype_len());
		fgets(line, STR_LEN_FILE_LINE, fp);
		h2.read(hd.allele_type().c_str(), line, hd.genotype_len());
		if (h1.length() != hd.genotype_len() || h2.length() != hd.genotype_len()) {
			Logger::error("Incorrect haplotype data for individual %d!", i);
			exit(1);
		}
		hd[i].setHaplotypes(h1, h2);
	}
	fclose(fp);
	hd.checkAlleleSymbol();
}

void HaploFile::writeHaploData(HaploData &hd, const char *suffix)
{
	FILE *fp;
	char output_file[STR_LEN_FILE_NAME], buf[STR_LEN_FILE_LINE], id;
	int i, j;
	m_haplo_data = &hd;
	makeFileName(output_file, suffix);
	fp = fopen(output_file, "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename);
		exit(1);
	}
	fprintf(fp, "%d\n", hd.genotype_num());
	fprintf(fp, "%d\n", hd.genotype_len());
	fprintf(fp, "P");
	for (i=0; i<hd.genotype_len(); i++) {
		fprintf(fp, " %d", hd.allele_postition(i));
	}
	fprintf(fp, "\n");
	fprintf(fp, "%s\n", hd.allele_type().c_str());
	for (i=0; i<hd.genotype_num(); i++) {
		id = hd[i].id()[0];
		if (id >= '0' && id <= '9') {
			fprintf(fp, "#%s\n", hd[i].id().c_str());
		}
		else {
			fprintf(fp, "%s\n", hd[i].id().c_str());
		}
		for (j=0; j<2; j++) {
			fprintf(fp, "%s\n", hd[i](j).write(hd.allele_type().c_str(), buf));
		}
	}
	fclose(fp);
}

void HaploFile::writePattern(HaploBuilder &hb, const char *suffix)
{
	FILE *fp;
	char output_file[STR_LEN_FILE_NAME], buf[STR_LEN_FILE_LINE];
	HaploPattern *hp;
	list<HaploPattern*>::iterator i_hp;
	m_haplo_data = hb.haplo_data();
	makeFileName(output_file, suffix);
	fp = fopen(output_file, "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", filename);
		exit(1);
	}
	fprintf(fp, "Frequency\tLength\t%s\n", writeAlleleName(buf));
	for (i_hp = hb.haplo_pattern().begin(); i_hp != hb.haplo_pattern().end(); i_hp++) {
		hp = *i_hp;
		fprintf(fp, "%f\t%d\t%s\n", hp->frequency() / m_haplo_data->genotype_num(), hp->length(), hp->write(buf, true));
	}
	fclose(fp);
}

char *HaploFile::readAlleleName(char *buffer)
{
	int i;
	char *s, *buf, *delim = " \t\r\n";
	buf = new char [strlen(buffer)+100];
	strcpy(buf, buffer);
	s = strtok(buf, delim);
	i = 0;
	while (s != NULL && i < m_haplo_data->genotype_len()) {
		m_haplo_data->setAlleleName(i, s);
		s = strtok(NULL, delim);
	}
	delete[] buf;
	return buffer;
}

char *HaploFile::writeAlleleName(char *buffer)
{
	int i;
	buffer[0] = 0;
	for (i=0; i<m_haplo_data->genotype_len(); i++) {
		strcat(buffer, m_haplo_data->allele_name(i).c_str());
		strcat(buffer, " ");
	}
	return buffer;
}

char *HaploFile::makeFileName(char *filename, const char *suffix)
{
	strcpy(filename, m_filename);
	if (suffix != NULL) {
		strcat(filename, suffix);
	}
	return filename;
}


////////////////////////////////
//
// class HaploFileHPM

void HaploFileHPM::readHaploData(HaploData &hd)
{
	FILE *fp;
	char line[STR_LEN_FILE_LINE];
	List<Haplotype> haplos;
	Haplotype *h, *h1, *h2;
	int i;
	m_haplo_data = &hd;
	fp = fopen(m_filename, "r");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename);
		exit(1);
	}
	fgets(line, STR_LEN_FILE_LINE, fp);
	checkHeader(line);
	while(fgets(line, STR_LEN_FILE_LINE, fp) != NULL) {
		h = new Haplotype;
		readHaplotype(*h, line);
		if (h->length() != hd.genotype_len()) {
			Logger::error("Incorrect haplotype data in line %d!", haplos.size()+2);
			exit(1);
		}
		haplos.addLast(h);
	}
	if (haplos.size() % 2 != 0) {
		Logger::error("Incorrect haplotype data in line %d!", haplos.size()+2);
		exit(1);
	}
	hd.setGenotypeNum(haplos.size()/2);
	haplos.moveFirst();
	for (i=0; i<hd.genotype_num(); i++) {
		h1 = haplos.getCurrent();
		haplos.moveNext();
		h2 = haplos.getCurrent();
		haplos.moveNext();
		hd[i].setHaplotypes(*h1, *h2);
		hd[i].setID(h1->id());
	}
	fclose(fp);
	hd.checkAlleleSymbol();
}

void HaploFileHPM::writeHaploData(HaploData &hd, const char *suffix)
{
	FILE *fp;
	char output_file[STR_LEN_FILE_NAME], buf[STR_LEN_FILE_LINE];
	int i, j;
	m_haplo_data = &hd;
	makeFileName(output_file, suffix);
	fp = fopen(output_file, "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename);
		exit(1);
	}
	fprintf(fp, "Id\t%s", writeAlleleName(buf));
	fprintf(fp, "\n");
	for (i=0; i<hd.genotype_num(); i++) {
		for (j=0; j<2; j++) {
			fprintf(fp, "%s\n", writeHaplotype(hd[i](j), buf));
		}
	}
	fclose(fp);
}

void HaploFileHPM::writeHaploDataWithFreq(HaploData &hd, const char *suffix)
{
	FILE *fp;
	char output_file[STR_LEN_FILE_NAME], buf[STR_LEN_FILE_LINE];
	int i, j;
	m_haplo_data = &hd;
	makeFileName(output_file, suffix);
	fp = fopen(output_file, "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename);
		exit(1);
	}
	fprintf(fp, "Id\t%s\tCONFIDENCE\n", writeAlleleName(buf));
	for (i=0; i<hd.genotype_num(); i++) {
		for (j=0; j<2; j++) {
			fprintf(fp, "%s\t%e\n", writeHaplotype(hd[i](j), buf), hd[i](j).weight());
		}
	}
	fclose(fp);
}

char *HaploFileHPM::readHaplotype(Haplotype &h, char *buffer)
{
	int i;
	char *s, *buf, *delim = " \t\r\n";
	buf = new char [strlen(buffer)+100];
	strcpy(buf, buffer);
	s = strtok(buf, delim);
	if (strlen(s) >= STR_LEN_GENOTYPE_ID) {
		Logger::warning("Too long Haplotype ID!");
		s[STR_LEN_GENOTYPE_ID-1] = 0;
	}
	h.setID(s);
	for (i=1; i<m_line_start; i++) {
		s = strtok(NULL, delim);
	}
	s += strlen(s)+1;
	s = h.read(NULL, s, m_haplo_data->genotype_len());
	if (m_weighted) {
		s += strspn(s, delim);
		h.setWeight(atof(s));
		s += strcspn(s, delim);
	}
	else {
		h.setWeight(1);
	}
	return s;
}

char *HaploFileHPM::writeHaplotype(const Haplotype &h, char *buffer)
{
	char *s = buffer;
	strcpy(s, h.id().c_str());
	strcat(s, "\t");
	s += strlen(s);
	h.write(NULL, s);
	return buffer;
}

void HaploFileHPM::checkHeader(char *buffer)
{
	int i, len;
	ConstStringList fields, names;
	const char *delim = " \t\r\n", *s;
	char *buf;
	fields.addLast("Id");
	fields.addLast("Status");
	fields.addLast("CONFIG_ID");
	fields.addLast("CONFIDENCE");
	buf = new char [strlen(buffer)+100];
	strcpy(buf, buffer);
	s = strtok(buf, delim);
	if (s != NULL && strcmp(s, "Id") != 0) {
		Logger::error("Invalid file type!");
		exit(1);
	}
	m_line_start = 1;
	s = strtok(NULL, delim);
	while (s != NULL && fields.contain(s)) {
		m_line_start++;
		s = strtok(NULL, delim);
	}
	len = 0;
	while (s != NULL && !fields.contain(s)) {
		len++;
		names.addLast(s);
		s = strtok(NULL, delim);
	}
	if (s != NULL && strcmp(s, "CONFIDENCE") == 0) {
		m_weighted = true;
	}
	else {
		m_weighted = false;
	}
	m_haplo_data->setGenotypeLen(len);
	s = names.getFirst();
	for (i=0; i<len; i++) {
		m_haplo_data->setAlleleName(i, s);
		s = names.getNext();
	}
	delete[] buf;
}


////////////////////////////////
//
// class HaploFileBench

HaploFileBench::HaploFileBench()
{
	m_children_file[0] = 0;
	m_posinfo_file[0] = 0;
	m_parents_num = 0;
	m_children_num = 0;
}

HaploFileBench::HaploFileBench(const char *filename, const char *posinfo, const char *children)
	: HaploFile(filename)
{
	if (children != NULL) strcpy(m_children_file, children);
	else m_children_file[0] = 0;
	if (posinfo != NULL) strcpy(m_posinfo_file, posinfo);
	else m_posinfo_file[0] = 0;
	m_parents_num = 0;
	m_children_num = 0;
}

void HaploFileBench::readHaploData(HaploData &hd)
{
	FILE *fp;
	char line[STR_LEN_FILE_LINE];
	char *s, *delim = " \t\r\n";
	List<Haplotype> haplos;
	Haplotype *h1, *h2;
	int i;
	m_haplo_data = &hd;
	fp = fopen(m_filename, "r");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename);
		exit(1);
	}
	// get genotype length
	fgets(line, STR_LEN_FILE_LINE, fp);
	s = line + strspn(line, delim);
	hd.setGenotypeLen(strcspn(s, delim));
	// set loci type
	for (i=0; i<hd.genotype_len(); i++) {
		hd.setAlleleType(i, 'S');
	}
	fclose(fp);
	readHaploFile(haplos, m_filename);
	m_parents_num = haplos.size() / 2;
	if (m_children_file[0] != 0) {
		readHaploFile(haplos, m_children_file);
		m_children_num = haplos.size() / 2 - m_parents_num;
	}
	hd.setGenotypeNum(m_parents_num + m_children_num);
	hd.setUnphasedNum(m_parents_num);
	haplos.moveFirst();
	for (i=0; i<hd.genotype_num(); i++) {
		h1 = haplos.getCurrent();
		haplos.moveNext();
		h2 = haplos.getCurrent();
		haplos.moveNext();
		hd[i].setHaplotypes(*h1, *h2);
		hd[i].setID(h1->id());
	}
	hd.checkAlleleSymbol();
}

void HaploFileBench::writeHaploData(HaploData &hd, const char *suffix)
{
	FILE *fp;
	char output_file[STR_LEN_FILE_NAME], buf[STR_LEN_FILE_LINE];
	int i, j;
	m_haplo_data = &hd;
	makeFileName(output_file, suffix);
	fp = fopen(output_file, "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename);
		exit(1);
	}
	for (i=0; i<hd.genotype_num(); i++) {
		for (j=0; j<2; j++) {
			fprintf(fp, "%s   %d 0 %s\n", writeHaplotype(hd[i](j), buf), i+2*j, hd[i](j).id().c_str());
		}
	}
	fclose(fp);
}

void HaploFileBench::writeHaploDataWithFreq(HaploData &hd, const char *suffix)
{
	FILE *fp;
	char output_file[STR_LEN_FILE_NAME], buf[STR_LEN_FILE_LINE];
	int i, j;
	m_haplo_data = &hd;
	makeFileName(output_file, suffix);
	fp = fopen(output_file, "w");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", m_filename);
		exit(1);
	}
	for (i=0; i<hd.genotype_num(); i++) {
		for (j=0; j<2; j++) {
			fprintf(fp, "%s   %f\n", writeHaplotype(hd[i](j), buf), hd[i](j).weight());
		}
	}
	fclose(fp);
}

void HaploFileBench::readHaploFile(List<Haplotype> &haplos, const char *filename)
{
	FILE *fp;
	char line[STR_LEN_FILE_LINE];
	char *s, *delim = " \t\r\n";
	Haplotype *h;
	int heterozygous;
	fp = fopen(filename, "r");
	if (fp == NULL) {
		Logger::error("Can not open file %s!", filename);
		exit(1);
	}
	// read haplotypes
	heterozygous = 1;
	while(fgets(line, STR_LEN_FILE_LINE, fp) != NULL) {
		h = new Haplotype;
		s = readHaplotype(*h, line, heterozygous);
		s += strspn(s, delim);		// next is number
		s += strcspn(s, delim);		// skip
		s += strspn(s, delim);		// next is 0
		s += strcspn(s, delim);		// skip
		s += strspn(s, delim);		// next is id
		s[strcspn(s, "\r\n")] = 0;
		h->setID(s);
		if (h->length() != m_haplo_data->genotype_len()) {
			Logger::error("Incorrect haplotype data in line %d of %s!", haplos.size()+2, filename);
			exit(1);
		}
		haplos.addLast(h);
		heterozygous = 3 - heterozygous;
	}
	if (haplos.size() % 2 != 0) {
		Logger::error("Incorrect haplotype data in line %d!", haplos.size()+2);
		exit(1);
	}
	fclose(fp);
}

char *HaploFileBench::readHaplotype(Haplotype &h, char *buffer, int heterozygous)
{
	int i;
	char *buf, *delim = " \t\r\n";
	h.setLength(m_haplo_data->genotype_len());
	buf = buffer + strspn(buffer, delim);
	for (i=0; i<h.length(); i++) {
		if (buf[i] == '0') {
			h[i] = -1;
		}
		else if (buf[i] == '9') {
			h[i] = '0' + heterozygous;
		}
		else {
			h[i] = (unsigned char) buf[i];
		}
	}
	buf += h.length();
	return buf;
}

char *HaploFileBench::writeHaplotype(const Haplotype &h, char *buffer)
{
	int i;

	for (i=0; i<h.length(); i++) {
		if (h[i].isMissing()) {
			buffer[i] = '0';
		}
		else {
			buffer[i] = h[i].asChar();
		}
	}
	buffer[h.length()] = 0;
	return buffer;
}
