#!/usr/bin/env ruby
# frozen_string_literal: true

require "nokogiri"
require "pathname"
require "uri"

site_dir = Pathname(ARGV.fetch(0, "_site")).expand_path
base_url = ARGV.fetch(1, "/AdditiveFOAM").sub(%r{/$}, "")
errors = []
documents = {}

site_dir.glob("**/*.html").each do |file|
  documents[file] = Nokogiri::HTML(file.read)
end

def output_file_for(site_dir, path)
  clean = path.sub(%r{^/}, "")
  candidate = site_dir.join(clean)
  return candidate if candidate.file?
  return candidate.join("index.html") if candidate.directory?
  return site_dir.join("#{clean}.html") if site_dir.join("#{clean}.html").file?

  candidate.join("index.html")
end

documents.each do |source_file, document|
  source_url = "/" + source_file.relative_path_from(site_dir).to_s
  source_url = source_url.sub(%r{/index\.html$}, "/")

  document.css("a[href]").each do |link|
    href = link["href"].to_s.strip
    next if href.empty? || href.match?(%r{^(https?:|mailto:|tel:|javascript:)})

    begin
      uri = URI.parse(href)
    rescue URI::InvalidURIError => e
      errors << "#{source_file}: invalid URI #{href.inspect}: #{e.message}"
      next
    end
    path = uri.path.to_s
    fragment = uri.fragment

    if path.empty?
      target_url = source_url
    elsif path.start_with?("/")
      unless path == base_url || path.start_with?("#{base_url}/")
        errors << "#{source_file}: root-relative link escapes baseurl: #{href}"
        next
      end
      target_url = path.delete_prefix(base_url)
      target_url = "/" if target_url.empty?
    else
      source_dir = source_url.end_with?("/") ? source_url : File.dirname(source_url) + "/"
      target_url = Pathname(source_dir).join(path).cleanpath.to_s
      target_url = "/#{target_url}" unless target_url.start_with?("/")
    end

    target_file = output_file_for(site_dir, target_url)
    unless target_file.file?
      errors << "#{source_file}: missing target #{href} (resolved to #{target_file})"
      next
    end

    next if fragment.nil? || fragment.empty?

    target_document = documents[target_file] || Nokogiri::HTML(target_file.read)
    anchor_found = target_document.css("[id], a[name]").any? do |node|
      node["id"] == fragment || node["name"] == fragment
    end
    errors << "#{source_file}: missing anchor #{href}" unless anchor_found
  end
end

if errors.empty?
  puts "Checked #{documents.length} HTML files: internal links are valid."
else
  warn errors.join("\n")
  exit 1
end
