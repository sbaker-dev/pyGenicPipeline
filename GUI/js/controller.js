
eel.expose(showStatus);
function showStatus(text, element_id, failed){
    /**
     * This will set element_id to have the text provided in green if successful and red if not
     *
     * @param {string} text: The message to be displayed
     * @param {string} element_id: The element to target
     * @param {boolean} failed: If this is a failure message or a success
     */

    const target_element = document.getElementById(element_id);
    target_element.innerHTML = text;
    target_element.classList.remove("text-danger", "text-success");
    if (failed){
        target_element.classList.add("text-danger")

    } else {
        target_element.classList.add("titles_foreground")
    }

}