#!/usr/bin/env python
# coding: utf-8

"""
send email
"""

# ---------
# Change Logs:
# 2016-12-22 Thursday 10:44    add multiple file attachment
# ---------

__author__ = 'Li Pidong'
__email__ = 'lipidong@126.com'
__version__ = '1.0.1'
__status__ = 'production'

import sys
reload(sys)
sys.setdefaultencoding('utf-8')
import argparse
import os

import smtplib
import email.mime.text
import email.header
from email.mime.application import MIMEApplication
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText


def send_email(to, header_text,
               message='This email was generated automatically.',
               cc=None, attachment=None, from_header='Panel', to_header=None):

    from_addr = 'lipidong@grandomics.com'
    password = 'Lipidong1990'
    smtp_server = 'smtp.exmail.qq.com'
    if cc:
        to_addr = to.split(',') + cc.split(',')
    else:
        to_addr = to.split(',')
    # msg = email.mime.text.MIMEText(message, 'text', 'utf-8')
    msg = MIMEMultipart()
    body = MIMEText(message, 'html', 'utf-8')
    msg.attach(body)
    msg['Subject'] = email.header.Header(header_text, 'utf-8')
    msg['From'] = email.header.Header(from_header, 'utf-8')
    if to_header:
        to_header = '; '.join(to_header.split(','))
        msg['To'] = email.header.Header(to_header, 'utf-8')
    else:
        msg['To'] = email.header.Header(to, 'utf-8')
    if cc:
        msg['Cc'] = cc
    if attachment:
        for xx in attachment.split(','):
            part = MIMEApplication(open(xx, 'rb').read())
            part.add_header('Content-Disposition', 'attachment',
                            filename=os.path.basename(xx))
            msg.attach(part)

    smtp = smtplib.SMTP()
    smtp.connect(smtp_server)
    smtp.login(from_addr, password)
    smtp.sendmail(from_addr, to_addr, msg.as_string())
    smtp.quit()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Send email')
    parser.add_argument('--to', help='recipients')
    parser.add_argument(
        '--header', help='the header of  email', default='autorun')
    parser.add_argument('--message', help='the body of email',
                        default='This email was generated automatically.')
    parser.add_argument(
        '--attachment', help='the attachment of email', default=None)
    parser.add_argument('--cc', help='cc', default=None)
    parser.add_argument(
        '--from_header', help='the from header of email', default='Bioinfo-auto')
    parser.add_argument(
        '--to_header', help='the to header of email', default=None)

    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    send_email(args.to, args.header, args.message, args.cc,
               args.attachment, args.from_header, args.to_header)
