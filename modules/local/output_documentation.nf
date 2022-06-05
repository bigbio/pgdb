/*
 * Output Description HTML
 */
process OUTPUT_DOCUMENTATION {

    container "nfcore/pgdb:1.0.0"
    
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs
    file images

    output:
    path 'results_description.html'

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}