/**
 * @file boundary/incident_field.h
 * @brief The incident field the scattered-field boundary adds: the TF/SF surface's E corrections
 *
 * The TF/SF surface separates the total-field from the scattered-field region; the E update of a
 * cell next to it adds the incident field's curl term. The boundary does not know where the incident
 * field comes from: the application hands it the source's term, sources::tfsf_e_terms, with
 * set_incident_e_terms, and cpml_tfsf_terms adds it after the CPML psi.
 */
#ifndef FDTD_BOUNDARY_INCIDENT_FIELD_H
#define FDTD_BOUNDARY_INCIDENT_FIELD_H

namespace fdtd::boundary
{
/// adds the incident-field term of component A at cell (i, j, k), scaled by the model's coefficient c
using IncidentETerms = void (*)(int A, int i, int j, int k, double c, double &e_new);
void set_incident_e_terms(IncidentETerms terms); ///< nullptr: no incident field
} // namespace fdtd::boundary

#endif
