#include <GkylVlasov.h>

Vlasov::Vlasov(unsigned cdim, unsigned vdim, unsigned polyOrder, const char* basisType, double qbym, bool hasForceTerm)
 : cdim(cdim), vdim(vdim), polyOrder(polyOrder), basisType(basisType), qbym(qbym), hasForceTerm(hasForceTerm) {
}

Vlasov::~Vlasov() {}

void Vlasov::setAuxFields(GkylCartField_t *emField) {
  emField = emField;
}

void* new_Vlasov(unsigned cdim, unsigned vdim, unsigned polyOrder, const char* basisType, double qbym, bool hasForceTerm) {
  Vlasov *vlasov = new Vlasov(cdim, vdim, polyOrder, basisType, qbym, hasForceTerm);
  return reinterpret_cast<void*>(vlasov);
}

void setAuxFields(Vlasov *eq, GkylCartField_t *emField) {
  eq->setAuxFields(emField);
}
