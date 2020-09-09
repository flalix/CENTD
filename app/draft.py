import logging
from logging.handlers import SMTPHandler, RotatingFileHandler
import os, sys, time, shutil

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
from app.dashPlots.layout    import layout
from app.dashPlots.callbacks import register_callbacks

from flask import Flask, request, current_app
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
from flask_login import LoginManager
from flask_mail import Mail
from flask_bootstrap import Bootstrap
from flask_moment import Moment
from flask_babel import Babel, lazy_gettext as _l
from flask.helpers import get_root_path
from flask_login import login_required
# from elasticsearch import Elasticsearch
from config import Config

sys.path.insert(1, '../src/')
from biobanco_lib   import *

db = SQLAlchemy()
migrate = Migrate()
login = LoginManager()
login.login_view = 'auth.login'
login.login_message = _l('Please log in to access this page.')
mail = Mail()
bootstrap = Bootstrap()
moment = Moment()
babel = Babel()
mydash = None

def create_app(config_class=Config):
    wantDash = True

    app = Flask(__name__)
    app.config.from_object(config_class)

    db.init_app(app)
    migrate.init_app(app, db)
    login.init_app(app)
    mail.init_app(app)
    bootstrap.init_app(app)
    moment.init_app(app)
    babel.init_app(app)

    app.file_index_css = rename_index_css()
    app.jinja_env.globals['file_index_css'] = app.file_index_css
    # print(">>> ", app.file_index_css)

    '''app.elasticsearch = Elasticsearch([app.config['ELASTICSEARCH_URL']]) \
        if app.config['ELASTICSEARCH_URL'] else None'''

    from app.errors import bp as errors_bp
    app.register_blueprint(errors_bp)

    from app.auth import bp as auth_bp
    app.register_blueprint(auth_bp, url_prefix='/auth')

    from app.main import bp as main_bp
    app.register_blueprint(main_bp)

    if not app.debug and not app.testing:
        if app.config['MAIL_SERVER']:
            auth = None
            if app.config['MAIL_USERNAME'] or app.config['MAIL_PASSWORD']:
                auth = (app.config['MAIL_USERNAME'],
                        app.config['MAIL_PASSWORD'])
            secure = None
            if app.config['MAIL_USE_TLS']:
                secure = ()
            mail_handler = SMTPHandler(
                mailhost=(app.config['MAIL_SERVER'], app.config['MAIL_PORT']),
                fromaddr='no-reply@' + app.config['MAIL_SERVER'],
                toaddrs=app.config['ADMINS'], subject='CENTD portal Failure',
                credentials=auth, secure=secure)
            mail_handler.setLevel(logging.ERROR)
            app.logger.addHandler(mail_handler)

        if app.config['LOG_TO_STDOUT']:
            stream_handler = logging.StreamHandler()
            stream_handler.setLevel(logging.INFO)
            app.logger.addHandler(stream_handler)
        else:
            if not os.path.exists('logs'):
                os.mkdir('logs')
            file_handler = RotatingFileHandler('logs/centd.log', maxBytes=10240, backupCount=10)
            file_handler.setFormatter(logging.Formatter(
                '%(asctime)s %(levelname)s: %(message)s '
                '[in %(pathname)s:%(lineno)d]'))
            file_handler.setLevel(logging.INFO)
            app.logger.addHandler(file_handler)

        app.logger.setLevel(logging.INFO)
        app.logger.info('CENTD portal startup')

    global mydash
    print("app __init__ mydash = ", mydash)
    mydash = init_dashboard_first(app)
    print("app __init__ mydash = ", mydash)

    '''
    global mydash
    mydash = init_dashboard_first(app)

    with app.app_context():
        # Import Dash application
        # from .plotlydash.dashboard import create_dashboard
        # dash_app = create_dashboard(app)
        global dboard
        dboard = Dashboard(app)
        return dboard.dash_app'''

    return mydash

def rename_index_css():
    want_rename = True

    if want_rename:
        _root = 'app/static/css/'
        filename = "index.css" # [x for x in os.listdir(_root) if 'index.css' in x][0]
        filenew = 'index_%s.css'%(time.ctime()).replace(' ','_').replace(":",'-')

        # print(">>> ", _root+filename, _root+filenew)
        shutil.copy(_root+filename, _root+filenew)
    else:
        filenew = 'index.css'

    return filenew


@babel.localeselector
def get_locale():
    return request.accept_languages.best_match(current_app.config['LANGUAGES'])


class Dashboard(object):
    def __init__(self, server):
        self.server   = server
        self.dash_app = self.create_dashboard(server)

    # https://hackersandslackers.com/plotly-dash-with-flask/
    def create_dashboard(self, server):
        """Create a Plotly Dash dashboard."""

        dash_app = dash.Dash(
            server=server,
            routes_pathname_prefix='/dashapp/',
            external_stylesheets=[
                '/static/css/styles.css',
            ]
        )

        # Create Dash Layout
        dash_app.layout = html.Div(id='dash-container')
        self.init_callbacks(dash_app)

        return dash_app

    def init_callbacks(self, dash_app):
        @dash_app.callback(Output('heatmap-id', 'figure'),
                          [Input('experiment-id', 'value')])
        def update_graph(vid):
            print("callback", vid)
            return vid

# https://medium.com/@olegkomarov_77860/how-to-embed-a-dash-app-into-an-existing-flask-app-ea05d7a2210b
def init_dashboard_first(app):

    # Meta tags for viewport responsiveness
    meta_viewport = {"name": "viewport", "content": "width=device-width, initial-scale=1, shrink-to-fit=no"}

    mydash = dash.Dash(__name__,
                         server=app,
                         url_base_pathname='/dashboard/',
                         assets_folder=get_root_path(__name__) + '/dashboard/assets/',
                         meta_tags=[meta_viewport])

    @mydash.callback(
        Output('heatmap-id', 'figure'),
        [Input('experiment-id', 'value')])
    def update_graph(vid):
        print(">>> update_graph", vid)
        return callback_update_graph(vid)

    with app.app_context():
        mydash.title = 'Dashapp 1'
        mydash.layout = layout
        register_callbacks(mydash)

    _protect_dashviews(mydash)

    return mydash

def callback_update_graph(vid, symmetric=True, want_nSamples=False):

    print(">>>> callback_update_graph", vid)
    vid = int(vid)

    try:
        from biobanco_lib import phm
    except:
        stri = "NOT FOUND - phm"
        print(stri)
        return None

    eperiment = phm.experiment
    dicexp    = phm.dicexp
    width     = phm.width

    height          = dicexp[vid]['height']
    name            = dicexp[vid]['name']
    fontsize        = dicexp[vid]['fontsize']
    sufix           = dicexp[vid]['author']
    filename        = dicexp[vid]['filename']
    experiment_type = dicexp[vid]['experiment_type']
    the_control     = dicexp[vid]['control']

    print(">>>> callback_update_graph eperiment:", eperiment, ' - ', name, "\n", filename)

    phm.bb = Biobanco(experiment_type, the_control, sufix, filename,
                      phm.root_data, phm.root_result, phm.exclude_exp,
                      nround=phm.nround, valpha=phm.valpha, verbose=phm.verbose)

    dfs, nSamples = phm.bb.create_heamap_table(verbose=False)

    if dfs is None:
        print("Nothing found")
        return None

    maxi = dfs.lfc.max()
    mini = dfs.lfc.min()

    if not symmetric:
        if maxi < _maxi: maxi = _maxi
        if mini > _mini: mini = _mini
    else:
        if abs(mini) > maxi:
            maxi = abs(mini)
        else:
            mini = -maxi

    template = 'plotly_white'; fontcolor='black'

    if want_nSamples:
        title = "LFC Heatmap: venoms x cytokines<br>%s model; cell type: '%s'<br>samples = %s"%(phm.bb.experiment_type, phm.bb.cell_type, nSamples)
    else:
        title = "LFC Heatmap: venoms x cytokines<br>%s model; cell type: '%s'"%(phm.bb.experiment_type, phm.bb.cell_type)

    return prepare_heatmap_fig(dfs, title, mini, maxi, width, height, template, fontsize, fontcolor)


def register_dashapps(app):
    from app.dashPlots.layout    import layout
    from app.dashPlots.callbacks import register_callbacks

    # Meta tags for viewport responsiveness
    meta_viewport = {"name": "viewport", "content": "width=device-width, initial-scale=1, shrink-to-fit=no"}

    global mydash
    mydash = dash.Dash(__name__,
                         server=app,
                         url_base_pathname='/dashboard/',
                         assets_folder=get_root_path(__name__) + '/dashboard/assets/',
                         meta_tags=[meta_viewport])

    with app.app_context():
        mydash.title = 'Dashapp 1'
        mydash.layout = layout
        register_callbacks(mydash)

    _protect_dashviews(mydash)

def _protect_dashviews(dashapp):
    for view_func in dashapp.server.view_functions:
        if view_func.startswith(dashapp.config.url_base_pathname):
            dashapp.server.view_functions[view_func] = login_required(dashapp.server.view_functions[view_func])


from app import models


def mutation_map():
    form = MutationMapForm()

    # dfh, countries = split_gisaid_fasta_headers(filehead, filefas, root_data, force=False)
    dfsumm = pdreadcsv(filesummary, root_protein, verbose=False)
    countries = list(dfsumm.country.unique())

    #--- metadata file
    dfmeta = pdreadcsv(filemeta, root_data0, verbose=False)


    #------------ objects ----------------------------------
    deschtml = Div(text=open("app/templates/mutmap.html").read(), sizing_mode="stretch_width")

    dicdata0 = {'country' : [], 'x' : [], 'y': [], 'aas': [], 'nns': []}
    dicdata1 = {'country' : [], 'x' : [], 'y': [], 'aas': [], 'nns': []}
    dicdata2 = {'country' : [], 'x' : [], 'y': [], 'aas': [], 'nns': []}

    source0 = ColumnDataSource(data=dicdata0)
    source1 = ColumnDataSource(data=dicdata1)
    source2 = ColumnDataSource(data=dicdata2)
    sources = [source0, source1, source2]

    increment = -1; ylabel = 'H (bits)'
    delwidth = .25

    colors = ['blue', 'red', 'green']

    key = 'S'
    countries2 = ['a', 'b', 'c' ]
    title0 = "Entropy for SARS-CoV-2 protein '%s'"
    title = title0%(key)

    TOOLTIPS=[
        ("country", "@country"),
        ("x", '@x{int}' ),
        ("aa", '@aas' ),
        ("#", "@nns"),
        ("h", '@y')
    ]

    #   min_border=0, toolbar_location=None, sizing_mode="scale_both"
    plot_width=1000; plot_height=600
    p = figure(title=title, plot_width=plot_width, plot_height=plot_height, tooltips=TOOLTIPS, sizing_mode="scale_both")

    # handles = []
    for i in range(3):
        handle = p.vbar(x='x', top='y', width=.25, color=colors[i], source=sources[i] ) # , legend=countries2[i]
        p.add_tools(HoverTool(renderers=[handle], tooltips=TOOLTIPS,  mode='vline')  )
        # handles.append(handle)


    p.x_range.range_padding = 0.6
    p.xgrid.grid_line_color = None
    # p.legend.location = "top_left"
    # p.legend.orientation = "vertical"

    xaxis = LinearAxis()
    # p.add_layout(xaxis, 'below') # major_label_overrides = xdic

    yaxis = LinearAxis()
    # p.add_layout(yaxis, 'left')

    p.add_layout(Grid(dimension=0, ticker=xaxis.ticker))
    p.add_layout(Grid(dimension=1, ticker=yaxis.ticker))

    p.xaxis.axis_label = xlabel
    p.yaxis.axis_label = ylabel

    # p.xaxis.major_label_overrides = xdic
    p.xaxis.major_label_orientation = math.pi/2

    # p.x_range=Range1d(0, xmaxi+5)
    # p.y_range=Range1d(0, maximus*1.2)

    p.title.text_font_size = '20pt'
    p.xaxis.axis_label_text_font_size = "14pt"
    p.yaxis.axis_label_text_font_size = "14pt"

    #---------------- controls ---------------------------------------
    # source_countries = ColumnDataSource(data={'country': countries})
    # def showmessage(title, content):
    #    dialogs(title, content, closable=True)

    txt_date_inf = TextInput(value="01/04/2020", title="Inferior (inclusive):")
    txt_date_sup = TextInput(value="12/05/2020", title="Superior (exclusive):")
    txt_dummy    = TextInput(value="---", title="dummy")
    txt_alert    = TextInput(value="---", title="alert")

    # min_month = Slider(title="initial month (2020)", start=1, end=6, value=4, step=1)
    # max_month = Slider(title="end month (2020)", start=1, end=7, value=6, step=1)

    '''
    checkbox_countries = CheckboxGroup(labels=countries, active=[0, 1, 2])
    '''
    autocomp0 = AutocompleteInput(completions=countries, value='Brazil',    background=colors[0])
    autocomp1 = AutocompleteInput(completions=countries, value='Argentina', background=colors[1])
    autocomp2 = AutocompleteInput(completions=countries, value='Chile',     background=colors[2])


    label_country0 = Div(text="""<br><b><span style="color: blue; ">Chose country:</span></b>""")
    label_country1 = Div(text="""    <b><span style="color: red;  ">Chose country:</span></b>""")
    label_country2 = Div(text="""    <b><span style="color: green;">Chose country:</span></b>""")
    group_countries = column(label_country0, autocomp0, label_country1, autocomp1, label_country2, autocomp2)

    '''
    reset_all = Button(label="reset all")

    def clear_checkbox_countries():
        checkbox_countries.active = [] * len(countries)

    reset_all.on_click(clear_checkbox_countries)

    group_countries = column(reset_all, checkbox_countries)
    '''
    keysall = keys + nspkeys
    # n = [i for i in range(len(keysall)) if keysall[i] = 'S'][0]
    select_protein = Select(options=keysall, value='S')

    label_protein = Div(text="""<br><b><span style="color: blue;">Chose a Protein/NSP</span></b><br>""")
    group_protein = column(label_protein, select_protein)

    label_date  = Div(text="""<b><span style="color: blue;">Chose a date interval</span></b><br>""")
    group_dates = column(label_date, txt_date_inf, txt_date_sup)

    button_submit = Button(label="refresh") # label='refresh',button_type='success'

    #---------------- udapte functions -------------------------------
    def update():

        _date_inf, _date_sup = select_ranges()
        if _date_inf is None: return

        print(_date_inf, _date_sup)

        global countries2
        '''
        countries2 = [checkbox_countries.labels[i] for i in checkbox_countries.active]
        if len(countries2) > 3:
            countries2 = countries2[:3]'''
        countries2 = [autocomp0.value, autocomp1.value, autocomp2.value]

        print(",".join(countries2))

        global key
        key = select_protein.value

        print("calc0", countries2[0])
        dicdata, nSamples0, maxVal0 = select_samples_by_country(countries2[0], 0, _date_inf, _date_sup)
        source0.data = dicdata

        print("calc1", countries2[1])
        dicdata, nSamples1, maxVal1 = select_samples_by_country(countries2[1], 1, _date_inf, _date_sup)
        source1.data = dicdata

        print("calc2", countries2[2])
        dicdata, nSamples2, maxVal2 = select_samples_by_country(countries2[2], 2, _date_inf, _date_sup)
        source2.data = dicdata

        maxVal = np.max( [maxVal0, maxVal1, maxVal2] )

        # .strftime('%d/%b/%Y')
        mat = [x + "(%d)" for x in countries2]
        stri = ", ".join(mat)%(nSamples0, nSamples1, nSamples2)

        p.title.text = "Protein '%s', countries: %s between %s and %s"%(key, stri, _date_inf, _date_sup )
        curdoc().title = "Protein '%s', countries: %s"%(key, ", ".join(countries2))


        '''
        try:
            del(hcirc0); del(hcirc1); del(hcirc2)
        except:
            pass
        '''

        '''
        hcirc0 = p.circle( x=0, y=0, radius=0, color=colors[0], legend_label= countries2[0])
        hcirc1 = p.circle( x=0, y=0, radius=0, color=colors[1], legend_label= countries2[1])
        hcirc2 = p.circle( x=0, y=0, radius=0, color=colors[2], legend_label= countries2[2])

        p.legend.location    = "top_left"
        p.legend.orientation = "vertical"
        '''

    '''
    def update_graph():
        p.legend.items = [ (countries2[0], handles[0]),
                           (countries2[1], handles[1]),
                           (countries2[2], handles[2])  ]

    '''

    def select_ranges():
        _date_inf = txt_date_inf.value
        _date_sup = txt_date_sup.value
        if len(_date_inf) != 10 or len(_date_sup) != 10:
            title = 'Format date'
            content = 'Date must be in the following format: dd/mm/yyyy'
            return None, None

        '''
        _date_inf = '01/%d/2020'%(min_month.value)
        _date_sup = '01/%d/2020'%(max_month.value)
        '''
        return _date_inf, _date_sup


    def select_samples_by_country(country, i, _date_inf, _date_sup, ini = 0, isProtein=True, increment = -1):
        keyName = key

        shan = Shannon(prjName, isProtein, country, keyName, base_filename,
                       root_data, root_protein, root_protent_figure, root_protent_entropy, root_protent_html, root_templates)

        #--- load metadata object
        dic = pdloaddic(shan.file_entropy, path=root_protent_entropy, verbose=False)
        # print(dic.keys())
        df_metadata = dic['metadata']

        fileprot = base_filename%(country, key)
        filefull = os.path.join(root_protein, fileprot)

        #--- load protein object
        mseqProt = MySequence(prjName, root_protein)
        mseqProt.readFasta(filefull, showmessage=False)

        #-- filter de i's related to the filter - convid_seqs_lib.py
        # print(df_metadata.shape, df_metadata.columns)
        ndxs = filter_hs_by_date_new(df_metadata, _date_inf, _date_sup, date_mask="%d/%m/%Y")
        # print(">>> ndxs", _date_inf, _date_sup, df_metadata.shape)
        # print(">>> ndxs", len(ndxs))
        mseqProt.seq_records = [mseqProt.seq_records[i] for i in ndxs]
        mseqProt.seqs        = [mseqProt.seqs[i]        for i in ndxs]
        nSamples  = len(mseqProt.seqs)
        print(">>> seqs", nSamples)
        #--- must recalculate all entropy again !
        ret = shan.fast_calcHSannon_bias_correction_1pos(mseqProt, df_metadata, type="Protein", verbose=False)

        matPiDic = shan.dicPiList
        matNDic  = shan.dicNList
        hs       = shan.HShannonList
        maxVal   = np.max(hs)

        if increment == -1:
            end = -1
        else:
            end = ini+increment

        maximus = []

        maxi = len(hs)
        if end == -1:
            end = maxi
        elif end > maxi:
            end = maxi

        matPiDic = matPiDic[ini : end]
        matNDic  = matNDic[ini : end]
        hs       = hs[ini : end]
        maximus.append(np.max(hs))

        labelx = []; _seqx = []; _seqn = []; count = 1
        for jj in range(len(matPiDic)):
            #---- percent
            dic2 = matPiDic[jj]
            val = "; ".join(["%s %.2f%%"%(k, dic2[k]*100) for k in dic2.keys() if dic2[k] != 0])
            _seqx.append(val)

            valks = []
            for k in dic2.keys():
                if dic2[k] == 0: continue

                if dic2[k] == 1:
                    valks = [k]
                    break

                valks.append("%s %.2f%%"%(k, dic2[k]*100) )

            labelx.append(str(count) + '-' + "; ".join(valks))

            #---- n or #
            dic3 = matNDic[jj]
            val = "; ".join(["%s %d"%(k, dic3[k]) for k in dic3.keys() if dic3[k] != 0])
            _seqn.append(val)

            count += 1

        seqx = np.arange(ini, end)
        if i > 0:
            seqx = [x + delwidth*i for x in seqx]

        xmaxi = len(hs)
        maximus = np.max(maximus)
        dicdata = {'country' : [country]*len(hs), 'x': seqx, 'y': hs, 'aas': _seqx, 'nns': _seqn}

        return dicdata, nSamples, maxVal


    def on_button_submit(event):
        if txt_dummy.value == 'a':
            txt_dummy.value = 'b'
        else:
            txt_dummy.value = 'a'

    button_submit.on_event(ButtonClick, on_button_submit)

    #---------------- ojbects ----------------------------------------
    # min_month, max_month,

    controls = [group_dates, button_submit, group_protein, group_countries]
    # txt_date_inf, txt_date_sup, checkbox_countries
    for control in [txt_dummy, txt_alert]:
        if control == txt_alert:
            # control.on_change('value', lambda attr, old, new: update())
            pass
        else:
            control.on_change('value', lambda attr, old, new: update())

    inputs = column(*controls, width=int(plot_width/4), height=plot_height)
    inputs.sizing_mode = "fixed"
    lout = layout([ [deschtml], [inputs, p] ], sizing_mode="scale_both" )

    update()  # initial load of the data
    # curdoc().add_root(lout)

    return render_template('mutation_map.html', title=_('Mutation Map - SARS-CoV-2'), form=form, lout=lout)

    
