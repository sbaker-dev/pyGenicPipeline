"""
This uses html templates and beautiful soup to construct the GUI. This is done programmatically to speed up
construction as each 'job' is basically the same page just with different arguments that belong to a different parent.
"""

from pyGenicPipeline import ArgMaker
import sys

try:
    from bs4 import BeautifulSoup
except ImportError:
    BeautifulSoup = None
    sys.exit("GUI Constructor requires bs4 BeautifulSoup")

from pathlib import Path


class GUIMaker(ArgMaker):
    def __init__(self):
        super().__init__()

    def construct_gui(self):

        # todo To compile, go to the eel page to get the commands to put into the command line

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
            title_base.find("div", attrs={"class": "container"}).extend(self._link_card(job_group, description))

        # Write the index controller
        # todo will be renamed to index when finished
        with open(Path(Path(__file__).parent, "Web", "TEST.html"), "w", encoding="utf-8") as file:
            file.write(str(title_base))

    def construct_job_type_page(self, job_group):
        # Construct the Job Type page
        base = open(Path(Path(__file__).parent, "Web", f"{job_group}.html"), "rb")
        job_type_base = BeautifulSoup(base, features="html.parser")

        # Create the Job Pages to Link to
        self._create_html_files(self.operations[job_group].keys())

        # Link Job Pages to the Job Type Page
        for job in self.operations[job_group]:
            if job != "Description":
                description = self.operations[job_group][job]["Description"]
                job_type_base.find("div", attrs={"class": "container"}).extend(self._link_card(job, description))
                self.construct_job_page(job, self.operations[job_group][job])

        with open(Path(Path(__file__).parent, "Web", f"{job_group}.html"), "w", encoding="utf-8") as file:
            file.write(str(job_type_base))

    def construct_job_page(self, job, job_args):

        base = open(Path(Path(__file__).parent, "Web", f"{job}.html"), "rb")
        job_base = BeautifulSoup(base, features="html.parser")

        for arg_type in job_args.keys():
            if arg_type != "Description":
                # todo Add separator elements to allow for the different types via parsing arg_type to the yaml file
                # todo Need to add a arg type seperator (this will allow us to generalise this)

                for arg in job_args[arg_type]:
                    description = self.yaml_parameters["Arg_Descriptions"][arg]
                    job_base.find("div", attrs={"class": "container"}).extend(self._arg_card(arg, description))

        with open(Path(Path(__file__).parent, "Web", f"{job}.html"), "w", encoding="utf-8") as file:
            file.write(str(job_base))

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

    def _link_card(self, job, info):
        """
        This updates a link card base with the information needed to link to the next job group / job.

        :param job: The current job name
        :param info: The job information, description, or other details to describe what is going to happen after
            clicking
        :return: BeautifulSoup data to extend the page this is to be placed on
        :rtype: BeautifulSoup
        """
        # Load cards
        card = open(Path(self.template_path, "link_card.html"), "rb")
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

    def _arg_card(self, arg, info):
        # Load cards
        card = open(Path(self.template_path, "arg_card.html"), "rb")
        card_base = BeautifulSoup(card, features="html.parser")

        # Set the Title ID and string based on the operation
        card_element = card_base.find("h2")
        card_element["id"] = arg
        card_element.string = arg.replace("_", " ")

        # Set the descriptions of the arg group
        description = card_base.find("p")
        description.string = info

        label = card_base.find("label")
        label["for"] = f"{arg}_input"

        text_input = card_base.find("input")
        text_input["id"] = f"{arg}_input"
        text_input["placeholder"] = f"Please enter the information for {arg.replace('_', ' ')} here"
        return card_base

    @property
    def template_path(self):
        """Root of templates"""
        return Path(Path(__file__).parent, "Web", "html_templates")


if __name__ == '__main__':
    GUIMaker().construct_gui()
