"""
This uses html templates and beautiful soup to construct the GUI. This is done programmatically to speed up
construction as each 'job' is basically the same page just with different arguments that belong to a different parent.
"""

from pyGenicPipeline import ArgMaker
from bs4 import BeautifulSoup
from pathlib import Path


class GUIMaker(ArgMaker):
    def __init__(self):
        super().__init__()

        print(self.yaml_parameters)

        print(self.operations)

    def construct_gui(self):

        # Construct the main page that we load into
        self.construct_title_page()

        # For each job type, construct a landing page to direct individuals to a specific job
        for job_type in self.operations:
            self.construct_job_type_page(job_type)

        return

    def construct_title_page(self):

        # Create the title page for this operation
        base = open(Path(self.template_path, "base.html"), "rb")
        title_base = BeautifulSoup(base, features="html.parser")

        self._create_html_files(self.operations.keys())

        # Add the job groupings that have been defined
        for job_group in self.operations:
            description = self.operations[job_group]["Description"]
            title_base.find("div", attrs={"class": "container"}).extend(self._linker_card(job_group, description))

        # Write the index controller
        # todo will be renamed to index when finished
        with open(Path(Path(__file__).parent, "Web", "TEST.html"), "w", encoding="utf-8") as file:
            file.write(str(title_base))

    def construct_job_type_page(self, job_group):
        # Base page
        base = open(Path(Path(__file__).parent, "Web", f"{job_group}.html"), "rb")
        job_base = BeautifulSoup(base, features="html.parser")

        self._create_html_files(self.operations[job_group].keys())

        for job in self.operations[job_group]:
            if job != "Description":
                description = self.operations[job_group][job]["Description"]
                job_base.find("div", attrs={"class": "container"}).extend(self._linker_card(job, description))

        with open(Path(Path(__file__).parent, "Web", f"{job_group}.html"), "w", encoding="utf-8") as file:
            file.write(str(job_base))

    def construct_job_page(self):
        return


    def _create_html_files(self, group_of_jobs):
        """
        Create the html pages we need for the pipeline dynamically

        :param group_of_jobs: A list of keys from a dict
        """

        # Load the base template
        base = open(Path(self.template_path, "base.html"), "rb")
        title_base = BeautifulSoup(base, features="html.parser")

        # For each page name within the group of jobs, removing descriptions of the current job, make another html page
        # as this job may have sub jobs vus have its own description
        for page in group_of_jobs:
            if page != "Description":
                with open(Path(Path(__file__).parent, "Web", f"{page}.html"), "w", encoding="utf-8") as file:
                    file.write(str(title_base))

    def _linker_card(self, job, info):
        """
        THis updates a linker card base with the information needed to link to the next job group / job.

        :param job: The current job name
        :param info: The job information, description, or other details to describe what is going to happen after
            clicking
        :return: BeautifulSoup data to extend the page this is to be placed on
        :rtype: BeautifulSoup
        """
        # Load cards
        card = open(Path(self.template_path, "card.html"), "rb")
        card_base = BeautifulSoup(card, features="html.parser")

        # Set the Title ID and string based on the operation
        card_element = card_base.find("h2")
        card_element["id"] = job
        card_element.string = job.replace("_", " ")

        # Set the descriptions of the job group
        description = card_base.find("p")
        description.string = info

        # Set the link to the next group of pages which will be a html named the job_group
        linker = card_base.find("a")
        linker["href"] = f"{job}.html"
        linker.string = f"Go to {job.replace('_', ' ')} operations"
        return card_base

    @property
    def template_path(self):
        """Root of templates"""
        return Path(Path(__file__).parent, "Web", "html_templates")


if __name__ == '__main__':
    GUIMaker().construct_gui()
