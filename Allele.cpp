
#include "Allele.h"

#include "MemLeak.h"


////////////////////////////////
//
// class AlleleSequence

bool AlleleSequence::isMatch(const AlleleSequence &as, int start1, int start2, int len) const
{
	int i;
	bool match;
	if (start1+len > length() || start2+len > as.length()) {
		match = false;
	}
	else {
		match = true;
		for (i=0; i<len; i++) {
			if (!m_alleles[start1+i].isMatch(as[start2+i])) {
				match = false;
				break;
			}
		}
	}
	return match;
}

int AlleleSequence::getDiffNum(const AlleleSequence &as, int start1, int start2, int len) const
{
	int i, match;
	match = 0;
	if (start1+len <= length() && start2+len <= as.length()) {
		for (i=0; i<len; i++) {
			if (!m_alleles[start1+i].isMatch(as[start2+i])) {
				match++;
			}
		}
	}
	return match;
}

int AlleleSequence::setLength(int len)
{
	if (len > 0) {
		m_alleles.resize(len);
	}
	else {
		m_alleles.clear();
	}
	return length();
}

char *AlleleSequence::readAllele(char type, char *buffer, Allele &allele)
{
	char *buf;
	buf = buffer;
	buf += strspn(buf, " \t\r\n");					// skip delimiters
	if (type == 'S') {							// SNP, bi-allelic
		if (buf[0] == '-' || buf[0] == '?') {	// missing allele
			allele = -1;
		}
		else {
			allele = (unsigned char) buf[0];
		}
		buf++;
	}
	else {										// microsatellite
		if (buf[0] == '-' || buf[0] == '?') {	// missing allele
			allele = -1;
		}
		else {
			allele = atoi(buf);
		}
		buf += strcspn(buf, " \t\r\n");
	}
	return buf;
}

char *AlleleSequence::writeAllele(char type, char *buffer, const Allele &allele)
{
	char buf[64];
	buffer[0] = 0;
	if (type == 'S') {							// SNP, bi-allelic
		if (allele.isMissing()) {				// missing allele
			buffer[0] = '?';
		}
		else {
			buffer[0] = allele.asChar();
		}
		buffer[1] = ' ';						// put delimiter
		buffer[2] = 0;
	}
	else {										// microsatellite
		sprintf(buf, "%d ", allele.asInt());
		strcat(buffer, buf);
	}
	return buffer;
}

char *AlleleSequence::read(const char *types, char *buffer, int len)
{
	int i;
	char *buf;
	Allele a;
	buf = buffer;
	if (types == NULL) {
		setLength(0);
		while (buf[0] != 0 && (len <=0 || length() < len)) {
			buf = readAllele('M', buf, a);
			(*this) += a;
		};
	}
	else {
		i = strlen(types);
		if (len > 0) {
			if (i < len) {
				Logger::error("Incomplete allele type information!");
				exit(1);
			}
		}
		else {
			len = i;
		}
		setLength(len);
		for (i=0; i<length(); ++i) {
			buf = readAllele(types[i], buf, m_alleles[i]);
		}
	}
	return buffer;
}

char *AlleleSequence::write(const char *types, char *buffer) const
{
	int i;
	char *buf;
	buffer[0] = 0;
	buf = buffer;
	if (types == NULL) {
		for (i=0; i<length(); i++) {
			writeAllele('M', buf, m_alleles[i]);
			buf = buffer + strlen(buffer);
		}
	}
	else {
		for (i=0; i<length(); i++) {
			writeAllele(types[i], buf, m_alleles[i]);
			buf = buffer + strlen(buffer);
		}
	}
	return buffer;
}
