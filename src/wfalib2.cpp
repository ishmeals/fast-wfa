#include "wfalib2.hpp"
// Include WFA2-lib header files
#include "external/WFA2-lib/wavefront/wavefront_align.h"

#include <string>

// WFA2-lib expects null-terminated C-strings, so we’ll convert std::string_view
int wfalib2_align(std::string_view a, std::string_view b, int x, int o, int e) {
    // Convert to std::string since wavefront_align needs a char*
    std::string pattern(a);
    std::string text(b);

    // Configure WFA alignment attributes for gap-affine mode
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine;
    attributes.affine_penalties.mismatch = x;
    attributes.affine_penalties.gap_opening = o;
    attributes.affine_penalties.gap_extension = e;
    attributes.alignment_scope = compute_scoret;

    // Initialize the WFA aligner
    wavefront_aligner_t* wf_aligner = wavefront_aligner_new(&attributes);

    // Align the sequences end-to-end
    wavefront_align(
        wf_aligner,
        pattern.c_str(), (int)pattern.size(),
        text.c_str(), (int)text.size()
    );

    // Retrieve the final alignment score
    int score = wf_aligner->cigar.score;

    // Clean up
    wavefront_aligner_delete(wf_aligner);

    return score;
}
