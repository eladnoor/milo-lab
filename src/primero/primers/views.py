# Create your views here.

from django.shortcuts import render_to_response
from django.http import HttpResponse, HttpResponseRedirect
from django.http import Http404
from models import Order, Primer, Pending
import time, sys, re
from forms import SearchForm, FetchMailForm
import email, imaplib, os
from pdfminer.pdfinterp import PDFResourceManager, process_pdf
from pdfminer.converter import TextConverter
from pdfminer.layout import LAParams
import StringIO
from forms import UploadSigmaPdfForm
from primero.primers.forms import UpdateLocationForm

def error(msg):
    return render_to_response('error.html', {'error_msg': msg})

def welcome(request):
    return HttpResponseRedirect('primero/search/')

def GetNowString():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime())
    
def primer_search(request):
    primers_list = Pending.objects.filter()
    return render_to_response('search.html', {'n_pending': len(primers_list)})

def fetch_login(request):
    return render_to_response('fetchlogin.html')

def example_error(request):
    return error("Example Error")

def results(request):
    search_form = SearchForm(request.GET)
    if not search_form.is_valid():
        primers_list = Primer.objects.all()
        query_string = ""
    else:
        query_string = search_form.cleaned_data["query"]
        try:
            primers_list = Primer.objects.filter(name__contains=query_string)
        except Primer.DoesNotExist:
            primers_list = []
    
    return render_to_response('results.html',
                              {'query': query_string,
                               'primers_list': primers_list})

def pending(request):
    try:
        primers_list = Pending.objects.all()
    except Primer.DoesNotExist:
        primers_list = []

    return render_to_response('pending.html',
                              {'primers_list': primers_list})

def update_location(request):
    loc_id = request.GET["loc_id"]
    query = request.GET["query"]
    primer = Primer.objects.get(loc_id=loc_id)

    if 'Delete' in request.POST:
        primer.delete()
        return HttpResponseRedirect("../query/?query=%s#" % query)
    elif 'Submit' in request.POST:
        update_form = UpdateLocationForm(request.POST)
        if not update_form.is_valid():
            return error("Update for is not valid")
        primer.box = update_form.cleaned_data["box"]
        if primer.box < 1 or primer.box > 99:
            return error("Box out of range: %d" % primer.box)
        primer.row = update_form.cleaned_data["row"]
        if ord(primer.row) < ord('A') or ord(primer.row) > ord('I'):
            return error("Row out of range: %s" % primer.row)
        primer.col = update_form.cleaned_data["col"]
        if primer.col < 1 or primer.col > 9:
            return error("Column out of range: %d" % primer.col)
        primer.loc_id = Primer.get_loc_id(primer.box, primer.row, primer.col)
        primer.save()
        #return render_to_response('updated.html',
        #    {'primers_list': primers_list, 'query': query})
        return HttpResponseRedirect("../query/?query=%s#" % query)
    else:    
        return render_to_response('change_location.html',
            {'primer': primer, 'query': query})

def update_pending(request):
    pending_list = []
    for check in request.POST:
        try:
            p = Pending.objects.filter(id=check)[0]
        except Exception: 
            continue
        pending_list.append(p)

    if 'Delete' in request.POST:
        for p in pending_list:
            p.delete()
        return HttpResponseRedirect('../pending/')
    elif 'Submit' in request.POST:
        updated_list = []
        for p in pending_list:
            existing_primers = Primer.objects.all()
            set_of_ids = set()
            for p in existing_primers:
                set_of_ids.add(p.loc_id)
            if len(set_of_ids) == 0:
                loc_id = 0
            else:
                all_ints_in_range = set(range(1, max(set_of_ids) + 2))
                available_ids = all_ints_in_range.difference(set_of_ids)
                loc_id = min(available_ids)
            
            # find plate, row, col
            box, row, col = Primer.loc_id2location(loc_id)
            entry = Primer(name=p.name, seq=p.seq,
                           comment='', box=box, row=row, col=col,
                           loc_id=loc_id)
            entry.save()
            p.delete()
            updated_list.append(entry)
        return render_to_response('updated.html', {'primers_list': updated_list})
    
    return error("Unknown Request")

def fetch_mail(request):
    fetch_form = FetchMailForm(request.POST)
    
    if not fetch_form.is_valid():
        return Http404
    
    user = fetch_form.cleaned_data["email"]
    pwd = fetch_form.cleaned_data["password"]
    detach_dir = '/tmp' # directory where to save attachments (default: current)
    
    # connecting to the gmail imap server
    m = imaplib.IMAP4_SSL("imap.gmail.com")
    resp, items = m.login(user, pwd)
    if resp != 'OK':
        return error("Wrong user name or password")
       
    #m.select("[Gmail]/INBOX") # here you a can choose a mail box like INBOX instead
    resp, items = m.select('sigma')
    if resp != 'OK':
        return error("Cannot locate label 'sigma'")
    
    resp, items = m.search(None, '(SUBJECT "Sigma-Aldrich Order confirmation")')
    items = items[0].split() # getting the mails id
    
    for emailid in items:
        # fetching the mail, "`(RFC822)`" means "get the whole stuff", but you can ask for headers only, etc
        resp, data = m.fetch(emailid, "(RFC822)")
        if data is None or data[0] is None:
            continue
        email_body = data[0][1] # getting the mail content
        mail = email.message_from_string(email_body) # parsing the mail content to get a mail object
    
        # Check if any attachments at all
        if mail.get_content_maintype() != 'multipart':
            continue
    
        print "["+mail["From"]+"] :" + mail["Subject"]
    
        # we use walk to create a generator so we can iterate on the parts and forget about the recursive headach
        for part in mail.walk():
            # multipart are just containers, so we skip them
            if part.get_content_maintype() == 'multipart':
                continue
    
            # is this part an attachment ?
            if part.get('Content-Disposition') is None:
                continue
    
            filename = part.get_filename()
            counter = 1
    
            # if there is no filename, we create one with a counter to avoid duplicates
            if not filename:
                filename = 'part-%03d%s' % (counter, 'bin')
                counter += 1
    
            att_path = os.path.join(detach_dir, filename)
    
            # Check if its already there
            if not os.path.isfile(att_path) :
                # finally write the stuff
                fp = open(att_path, 'wb')
                fp.write(part.get_payload(decode=True))
                fp.close()
                apply_lbl_msg = m.uid('COPY', emailid, 'done')
                if apply_lbl_msg[0] == 'OK':
                    mov, data = m.uid('STORE', emailid , '+FLAGS', '(\Deleted)')
                    m.expunge()                
                for primer in _pdfMine(att_path):
                    entry = Pending(name=primer[0], seq=primer[1])
                    entry.save()
    return HttpResponse("Done")

def upload_sigma_pdf(request):
    if request.method == 'POST':
        form = UploadSigmaPdfForm(request.POST, request.FILES)
        if form.is_valid():
            primers_list = []
            for primer in pdfMine(request.FILES['file']):
                entry = Pending(name=primer[0], seq=primer[1])
                entry.save()
                primers_list.append(entry)
            return render_to_response('upload_success.html',
                                      {'primers_list': primers_list})
        else:
            return error("Invalid operation")
    else:
        form = UploadSigmaPdfForm()
    return render_to_response('upload.html', {'form': form})

def _pdfMine(fname):
    fp = file(fname, 'rb')
    pdfMine(fp, sys.stdout)
    fp.close()

def pdfMine(fp):
    str_io = StringIO.StringIO()
    laparams = LAParams()
    
    rsrcmgr = PDFResourceManager(caching=True)
    device = TextConverter(rsrcmgr, str_io, codec='utf-8', laparams=laparams)
    
    process_pdf(rsrcmgr, device, fp, pagenos=set(), maxpages=0, password='',
                caching=True, check_extractable=True)
    
    pdf_string = str_io.getvalue()
    blocks = pdf_string.split('VC00')
    primers = []
    for b in blocks:
        for x in re.findall('EA[\s\n\d\.\W]+(\D+[\w\d\s\S]+)\nUMO', b):
            primers.append(x.split('\n'))
    return primers
