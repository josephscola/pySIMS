def _str2float(s: str) -> float:
    return float(s or 0)


class Semantic:
    def __init__(self):
        pass

    def meta_section_line(self, ast):
        """
        Semantic for a line of a section of metadata.
        
        :return: format the parsed line as [param, value]
        :rtype: list
        """
        return [ast.param, ast.value]

    def meta_subsection(self, ast):
        """
        Semantic for a subsection of metadata

        :return: a list with the name of the subsection and a dict of
                 parameters and their values
        :rtype: tuple
        """
        return [ast.name, dict(ast.lines)]

    def meta_section(self, ast):
        """
        Semantic for a section of metadata

        :return: a dictionnary containing the parameters with their
                 values and the different subsection : the subsection's name
                 is a key and their value is a dict of their
                 parameters/values
        :rtype: dict
        """
        params = dict(ast.body.lines)
        if ast.body.subsections:
            subsections = dict(ast.body.subsections)
            params.update(subsections)
        return params

    def calib_species_subsection(self, ast):
        """
        Semantic for the 'species' subsections in the "CALIBRATION
        PARAMETER" section.  The subsection's name is the name of the
        species

        :return: a list containing the name of the species and a dict
                 containing the associated parameters
        :rtype: list
        """
        return [ast.name, dict(ast.params)]

    def calib_param_section(self, ast):
        """
        Semantic for the section "CALIBRATION PARAMETER".

        :return: a dictionnary of parameters
        :rtype: dict
        """
        section_name = ast.header.lower()
        params = dict(ast.body.lines)
        species = {section_name: dict(ast.body.species)}
        params.update(species)
        return params

    def meta_csv_section(self, ast):
        """
        Semantic for section of metadata which are in a csv format

        :return: a dictionnary containing the name of the section and
                 the dv data as a 2D array.
        :rtype: dict
        """
        section_name = ast.header.section_name.lower()
        data = ast.data
        return {section_name: data}

    def data_section(self, ast):
        """
        Semantic for the main data section

        :return: a dict of elements.  Each elements is a dict
                 containing lists of data (floats).
        :rtype: dict
        """
        body = ast.body
        raw_data = [list(map(_str2float, row)) for row in body.data]

        data_dict = dict()
        data_header = body.data_header
        for i, e in enumerate(body.table_header):
            elem = e[0]
            nb_columns = len(e[1])
            values = lambda j: [row[i*nb_columns + j] for row in raw_data]
            data_dict[elem] = {data_header[i*nb_columns + j]: values(j) for j in range(nb_columns)}
        return data_dict

    def start(self, ast):
        """
        Semantic for the "start" section (the entry point of the parser).  It
        is used to aggregates all the different data and medata parsed
        from the different sections.

        :return: a dict containing the data and metadata
        :rtype: dict
        """
        sections = ast.sections
        metadata = {}
        data = {}
        for section in sections:
            if "metadata" in section:
                metadata.update(section["metadata"])
            if "data" in section:
                data.update(section["data"])
        return {"data": data, "metadata": metadata}
